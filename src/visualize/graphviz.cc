#include "graphviz.h"
using namespace std;

DEFINE_string(graphviz_engine, "fdp", "fdp (default), dot, neato, twopi, circo, sfdp, ...");

namespace agl {
//
// Input
//

graphviz::graphviz(const G &g) : graphviz_engine_(FLAGS_graphviz_engine) {
  for (V v : g.vertices()) {
    vertex_property_[v] = property_map();
  }

  for (V v : g.vertices()) {
    for (E e : g.edges(v)) {
      edge_property_[make_pair(v, to(e))]["dir"] = "forward";
    }
  }
}

void graphviz::read_dot(const string &filename) {
  ifstream ifs(filename);
  CHECK_PERROR(ifs);
  read_dot(ifs);
}

void graphviz::read_dot(istream &is) {
  // TODO: more robust implementation
  clear();

  // Header
  {
    string s;
    CHECK(is >> s && (s == "graph" || s == "digraph"));
    CHECK(is >> s);
    if (s != "{") {
      CHECK(is >> s && s == "{");
    }
  }
  // Commands
  while (!is.eof()) {
    string cmd = "";
    // TODO: ';' in escaped strings
    for (char c; !is.eof() && is.get(c) && c != ';'; ) cmd += c;

    vector<string> tokens = parse_space_separated_string<string>(cmd);
    if (tokens.size() < 2) continue;

    // "hoge=piyo,foo=bar"
    string property_str;
    {
      auto v = split(cmd, '[');
      CHECK(v.size() == 2);
      property_str = v[1];
      CHECK(property_str.back() == ']');
      property_str.pop_back();
    }

    // Where to put the properties
    property_map *target_property_map = nullptr;
    {
      if (tokens[1] == "->") {
        target_property_map = &edge_property_[make_pair(stoi(tokens[0]), stoi(tokens[2]))];
      } else {
        if (tokens[0] == "graph") {
          target_property_map = &graph_property_;
        } else if (tokens[0] == "node") {
          target_property_map = &vertex_common_property_;
        } else {
          target_property_map = &vertex_property_[stoi(tokens[0])];
        }
      }
    }

    // [name1, value1, name2, value2, ...]
    vector<string> name_and_values;
    {
      // TODO: '~' in escaped strings
      auto vs = split(property_str, '=');
      for (auto j : make_irange(vs.size())) {
        const string &s = vs[j];
        string::size_type i = s.rfind(',');
        if (j + 1 == vs.size() || i == string::npos) {
          name_and_values.emplace_back(strip(s));
        } else {
          name_and_values.emplace_back(strip(s.substr(0, i)));
          name_and_values.emplace_back(strip(s.substr(i + 1)));
        }
      }
    }

    CHECK(name_and_values.size() % 2 == 0);
    for (size_t i : make_irange(name_and_values.size() / 2)) {
      string &name = name_and_values[i * 2];
      string &value = name_and_values[i * 2 + 1];
      if (value[0] == '"') {
        CHECK(value.length() >= 2 && value.back() == '"');
        value = value.substr(1, value.length() - 2);
      }
      target_property_map->insert(make_pair(name, value));
    }
  }
}

void graphviz::clear() {
  graphviz_engine_ = FLAGS_graphviz_engine;
  graph_property_.clear();
  vertex_common_property_.clear();
  vertex_property_.clear();
  edge_property_.clear();
}

//
// Output
//
void graphviz::output_dot(ostream &os) const {
  os << "digraph G {" << endl;

  os << "  graph ";
  output_property_string(os, graph_property_);
  os << ";" << endl;

  for (const auto &i : vertex_property_) {
    os << "  " << i.first << " ";
    output_property_string(os, i.second);
    os << ";" << endl;
  }

  for (const auto &i : edge_property_) {
    os << "  " << i.first.first << " -> " << i.first.second << " ";
    output_property_string(os, i.second);
    os << ";" << endl;
  }

  os << "}" << endl;
}

void graphviz::output_dot(const string &filename) const {
  ofstream ofs(filename.c_str());
  CHECK_PERROR(ofs);
  output_dot(ofs);
  CHECK_PERROR(ofs);
}

void graphviz::generate_image(const string &img_filename, const string &img_type) const {
  string dot_filename = img_filename + ".dot";
  output_dot(dot_filename);

  string cmd = graphviz_engine_ + string(" ") + graphviz_option_ + " " + dot_filename + " -T " + img_type + " -o " + img_filename;
  CHECK(0 == system(cmd.c_str()));
}

void graphviz::generate_png(const string &filename) const {
  generate_image(filename, "png");
}

void graphviz::generate_eps(const string &filename) const {
  generate_image(filename, "eps");
}


//
// Direct property manipulation
//

void graphviz::set_graph_property(const std::string& property,
                                  const std::string& value) {
  graph_property_[property] = value;
}

void graphviz::set_vertex_property(V v, const std::string& property,
                                   const std::string value) {
  vertex_property_[v][property] = value;
}

void graphviz::set_vertex_property(
    const std::string& property,
    const std::unordered_map<V, std::string>& value_map) {
  set_vertex_property(property, [&](V v) -> string { return value_map.at(v); });
}

void graphviz::set_vertex_property(const std::string& property,
                                   std::function<std::string(V)> value_func) {
  for (auto &i : vertex_property_) {
    i.second[property] = value_func(i.first);
  }
}

void graphviz::set_edge_property(V u, V v, const std::string& property,
                                 const std::string value) {
  if (edge_property_.count({u, v}) == 0) swap(u, v);
  CHECK(edge_property_.count({u, v}));
  edge_property_[{u, v}][property] = value;
}

void graphviz::set_edge_property(
    const std::string& property,
    const std::unordered_map<std::pair<V, V>, std::string>& value_map) {
  set_edge_property(property, [&](V u, V v) -> string { return value_map.at({u, v}); });
}

void graphviz::set_edge_property(const std::string& property,
                                 std::function<std::string(V, V)> value_func) {
  for (auto &i : edge_property_) {
    i.second[property] = value_func(i.first.first, i.first.second);
  }
}

//
// Other configuration
//
void graphviz::merge_reverse_edges() {
  vector<pair<V, V>> edges_to_be_removed;

  for (const auto &i : edge_property_) {
    V v0 = i.first.first, v1 = i.first.second;
    if (edge_property_.count({v1, v0}) &&
        make_pair(v0, v1) > make_pair(v1, v0)) {
      edges_to_be_removed.emplace_back(v0, v1);
      edge_property_[{v1, v0}]["dir"] = "both";
    }
  }

  for (auto e : edges_to_be_removed) {
    edge_property_.erase(e);
  }
}

void graphviz::ignore_direction() {
  merge_reverse_edges();
  set_edge_property("dir", [](V, V) { return "none"; });
}

void graphviz::ignore_isolated_vertex(const G &g) {
  for (V v : g.vertices()) {
    if (g.degree(v, kFwd) + g.degree(v, kBwd) == 0) vertex_property_.erase(v);
  }
}

//
// Utility functions
//
void graphviz::output_property_string(std::ostream &os, const property_map &pm) {
  os << "[";
  bool f = true;
  for (const auto &p : pm) {
    if (f) f = false;
    else os << ", ";
    os << p.first << " = \"" << p.second << "\"";
  }
  os << "]";
}


//
// Easy functions
//
void graphviz_draw_graph(const G &g, const char *filename, const char *engine) {
  graphviz gz(g);
  gz.set_graphviz_engine(engine);

  gz.set_graph_property("splines", "true");
  gz.set_graph_property("overlap", "scale");
  // gz.set_vertex_property("label", [](V) { return ""; });
  gz.set_vertex_property("shape", [](V) { return "circle"; });

  for (V u : g.vertices()) for (V v : g.neighbors(u)) gz.set_edge_property(u, v, "label", "");

  gz.ignore_isolated_vertex(g);
  gz.ignore_direction();


  gz.generate_png(filename);
}

void graphviz_draw_edge_centrality(const G& g, const edge_centrality_map& ec,
                                   const char* png_filename,
                                   const char* graphviz_engine) {
  double min_value = numeric_limits<double>::max();
  double max_value = numeric_limits<double>::min();
  for (const auto &i : ec) {
    min_value = min(min_value, i.second);
    max_value = max(max_value, i.second);
  }

  graphviz gz(g);

  gz.set_graph_property("splines", "true");
  gz.set_graph_property("overlap", "scale");
  gz.set_vertex_property("shape", [](V) { return "point"; });
  gz.set_vertex_property("label", [](V) { return ""; });

  for (V u : g.vertices()) {
    for (V v : g.neighbors(u)) {
      const double x = ec.at({u, v});
      // gz.set_edge_property(u, v, "label", to_string(x));

      double A = 1.0 / 3.0;
      char buf[256];
      sprintf(buf, "%.3f 1.000 1.000", 1.0 - (((x - min_value) / (max_value - min_value)) * (1.0 - A) + A));
      gz.set_edge_property(u, v, "color", buf);
    }
  }

  gz.ignore_isolated_vertex(g);
  gz.ignore_direction();
  gz.set_graphviz_engine(graphviz_engine);
  gz.generate_png(png_filename);
}
}  // namespace agl
