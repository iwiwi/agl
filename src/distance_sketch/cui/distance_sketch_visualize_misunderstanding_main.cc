#include "easy_cui.h"
#include "distance_sketch/distance_sketch.h"
using namespace distance_sketch;

DEFINE_string(dir, "vis", "");

namespace {
int image_id = 0;

string double_to_string(double x) {
  ostringstream oss;
  oss << fixed << setprecision(2) << x;
  return oss.str();
}

string int_to_string_02d(int i) {
  ostringstream oss;
  oss << setfill('0') << setw(2) << right << i;
  return oss.str();
}

void generate_original_graph_layout(const G &g, const rank_array &ranks) {
  graphviz gz(g);
  gz.set_graph_property("splines", "true");
  gz.set_graph_property("overlap", "scale");
  gz.set_graph_property("size", "20");
  gz.set_graph_property("rankdir", "LR");
  gz.ignore_isolated_vertex(g);

  gz.set_vertex_property("shape", [](V) { return "Mrecord"; });
  gz.set_vertex_property("label", [&](V v) {
    return to_string(v + 1) + "|" + double_to_string(ranks[v] / (double)numeric_limits<rank_type>::max());
    // return "{" + to_string(v) + "|" + double_to_string(ranks[v] / (double)numeric_limits<rank_type>::max()) + "}";
  });

  // gz.ignore_direction();
  gz.generate_image(FLAGS_dir + "/layout.dot", "dot");
}

struct distance_sketch_visualizer {
  graphviz gz;

  distance_sketch_visualizer() {
    gz.read_dot(FLAGS_dir + "/layout.dot");
    gz.set_graphviz_option("-n");
    gz.set_graphviz_engine("neato");
  }

  void add_sketch_edges(const all_distances_sketches *ads, string color, bool draw_weight = false) {
    for (V v : make_irange(ads->sketches.size())) {
      for (auto e : ads->sketches[v]) {
        if (e.d > 0) {
          auto p = make_pair(min(v, e.v), max(v, e.v));
          gz.add_edge(p.first, p.second);

          string dir_str = gz.edge_properties()[p]["dir"];
          int dir = dir_str == "both" ? 3 :
              dir_str == "forward" ? 1 :
              dir_str == "back" ? 2 :
              0;
          dir |= (p.first == v ? 1 : 2);
          dir_str = dir == 3 ? "both":
              dir == 2 ? "back" :
              dir == 1 ? "forward":
              "none";

          gz.set_edge_property(v, e.v, "dir", dir_str);
          gz.set_edge_property(v, e.v, "color", color);
          gz.set_edge_property(v, e.v, "style", "solid");
          if (e.d > 1) gz.set_edge_property(v, e.v, "label", to_string(e.d));
        }
      }
    }
  }

  void generate() {
    for (const char *image_type : {"png", "eps"}) {
      gz.generate_image(FLAGS_dir + "/" + int_to_string_02d(image_id) + "." + image_type,
                        image_type);
    }
    ++image_id;
  }

  void highlight_vertices_in_sketch(const vertex_sketch_raw &s,
                          const string &property, const string &value) {
    for (const auto &e : s) {
      gz.set_vertex_property(e.v, property, value);
    }
  }
};
}  // namespace

int main(int argc, char **argv) {
  CHECK(0 == system(("mkdir -p " + FLAGS_dir).c_str()));

  JLOG_INIT(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);

  G g({
      {0, 1}, {1, 2}, {0, 3}, {3, 4}, {4, 5},
    });

  // G g(generate_grid(1, 4)); // easy_cui_init(argc, argv);
  vector<pair<V, V>> original_es = g.edge_list();
  sort(original_es.begin(), original_es.end());

  vector<double> my_ranks = {
    0.93, 0.60, 0.24, 0.49, 0.37, 0.02
  };

  rank_array ranks(my_ranks.size());
  for (size_t i = 0; i < my_ranks.size(); ++i) {
    ranks[i] = my_ranks[i] * numeric_limits<unsigned long long>::max();
  }

  auto ads = compute_all_distances_sketches(g, FLAGS_distance_sketch_k, ranks);
  auto srs = compute_sketch_retrieval_shortcuts(g, FLAGS_distance_sketch_k, ranks);

  generate_original_graph_layout(g, ranks);


  for (int b : make_irange(4)) {
    distance_sketch_visualizer dsv;
    if (b) dsv.gz.set_edge_property("style", [](V, V){ return "dashed"; });
    if (b & 1) dsv.add_sketch_edges(&ads, "orange");
    if (b & 2) dsv.add_sketch_edges(&srs, "blue", true);
    dsv.gz.set_edge_property("arrowsize", [](V, V) { return "0.5"; });
    dsv.generate();
  }

  for (V v_source : g.vertices()) {
    distance_sketch_visualizer dsv;
    const auto &s = ads.sketches[v_source];
    dsv.highlight_vertices_in_sketch(s, "style", "filled");
    dsv.highlight_vertices_in_sketch(s, "fillcolor", "gold");
    dsv.gz.set_vertex_property(v_source, "fillcolor", "pink");
    dsv.gz.set_edge_property("arrowsize", [](V, V) { return "0.7"; });
    dsv.generate();

    dsv.gz.set_edge_property("style", [](V, V){ return "dashed"; });
    dsv.add_sketch_edges(&srs, "lightskyblue");

    vector<int> ds(g.num_vertices(), -1000);
    for (auto e : s) ds[e.v] = e.d;

    dsv.gz.set_edge_property("style", [](V, V) { return "dashed"; });

    for (auto &i : dsv.gz.edge_properties()) {
      V u = i.first.first, v = i.first.second;
      if (abs(ds[v] - ds[u]) == 1) {
          i.second["penwidth"] = "2";
          i.second["style"] = "solid";
          i.second["color"] = "blue";
        }
      if (ds[v] == ds[u] + 1) i.second["dir"] = "forward";
      if (ds[u] == ds[v] + 1) i.second["dir"] = "back";
    }
    for (V v : g.vertices()) {
      for (auto e : srs.sketches[v]) {
        if (ds[e.v] == ds[v] + e.d) {
          dsv.gz.set_edge_property(v, e.v, "penwidth", "2");
          dsv.gz.set_edge_property(v, e.v, "style", "solid");
          dsv.gz.set_edge_property(v, e.v, "color", "blue");
        }
      }
    }

    dsv.gz.set_edge_property("arrowsize", [](V, V) { return "0.7"; });

    dsv.generate();
  }

  return 0;
}
