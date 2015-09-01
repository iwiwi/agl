#include "graphviz.h"
using namespace std;

DEFINE_string(graphviz_engine, "fdp", "fdp (default), dot, neato, twopi, circo, sfdp, ...");

namespace agl {
graphviz::graphviz(const G &g) :
    g_(g), graphviz_engine_(FLAGS_graphviz_engine) {
  for (V v : g_.vertices()) {
    vertex_property_[v] = property_map();
  }

  for (V v : g_.vertices()) {
    for (E e : g_.edges(v)) {
      edge_property_[make_pair(v, to(e))]["dir"] = "forward";
    }
  }
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

void graphviz::generate_iamge(const string &img_filename, const string &img_type) const {
  string dot_filename = img_filename + ".dot";
  output_dot(dot_filename);

  string cmd = graphviz_engine_ + string(" ") + dot_filename + " -T " + img_type + " -o " + img_filename;
  CHECK(0 == system(cmd.c_str()));
}

void graphviz::generate_png(const string &filename) const {
  generate_iamge(filename, "png");
}

void graphviz::generate_eps(const string &filename) const {
  generate_iamge(filename, "eps");
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

void graphviz::ignore_isolated_vertex() {
  for (V v : g_.vertices()) {
    if (g_.degree(v) == 0) vertex_property_.erase(v);
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

  gz.ignore_isolated_vertex();
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
  gz.set_vertex_property("shape", [](V) { return "circle"; });

  for (V u : g.vertices()) {
    for (V v : g.neighbors(u)) {
      const double x = ec.at({u, v});
      gz.set_edge_property(u, v, "label", to_string(x));

      double A = 1.0 / 3.0;
      char buf[256];
      sprintf(buf, "%.3f 1.000 1.000", 1.0 - (((x - min_value) / (max_value - min_value)) * (1.0 - A) + A));
      gz.set_edge_property(u, v, "color", buf);
    }
  }

  gz.ignore_isolated_vertex();
  gz.ignore_direction();
  gz.set_graphviz_engine(graphviz_engine);
  gz.generate_png(png_filename);
}
}  // namespace agl
