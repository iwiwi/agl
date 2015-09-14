#include "easy_cui.h"
#include "distance_sketch/distance_sketch.h"
using namespace distance_sketch;

namespace {
string double_to_string(double x) {
  ostringstream oss;
  oss << fixed << setprecision(2) << x;
  return oss.str();
}
}  // namespace

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);
  vector<pair<V, V>> original_es = g.edge_list();
  sort(original_es.begin(), original_es.end());

  rank_array ranks = generate_rank_array(g.num_vertices());
  auto ads = compute_all_distances_sketches(g, FLAGS_distance_sketch_k, ranks);
  auto srs = compute_sketch_retrieval_shortcuts(g, FLAGS_distance_sketch_k, ranks);


  {
    graphviz gz(g);
    gz.set_graph_property("splines", "true");
    gz.set_graph_property("overlap", "scale");
    gz.set_graph_property("size", "20");
    gz.ignore_isolated_vertex(g);

    gz.set_vertex_property("shape", [](V) { return "Mrecord"; });
    gz.set_vertex_property("label", [&](V v) {
      return "{" + to_string(v) + "|" + double_to_string(ads.ranks[v] / (double)numeric_limits<rank_type>::max()) + "}";
    });

    gz.ignore_direction();
    gz.generate_image("tmp.dot", "dot");
    gz.generate_png("01_graph.png");
  }

  graphviz gz;
  gz.read_dot("tmp.dot");
  //gz.set_graphviz_option("-s1 -n100");
  gz.set_graphviz_option("-n");

  map<pair<V, V>, int> map_ads, map_srs;
  for (V v : g.vertices()) {
    {
      for (auto e : ads.sketches[v]) {
        if (e.d > 1) gz.add_edge(v, (int)e.v);
        map_ads.emplace(make_pair(v, (int)e.v), (int)e.d);
      }
    }
    {
      for (auto e : srs.sketches[v]) {
        if (e.d > 1) gz.add_edge(v, (int)e.v);
        map_srs.emplace(make_pair(v, (int)e.v), (int)e.d);
      }
    }
  }
  gz.set_edge_property("color", [&](V u, V v) {
    return binary_search(original_es.begin(), original_es.end(), make_pair(u, v)) ? "black" :
        map_srs.count(make_pair(u, v)) ? "blue" : "red";
  });

  for (auto &i : gz.edge_properties()) {
    if (binary_search(original_es.begin(), original_es.end(), i.first)) {
      i.second["dir"] = "none";
    }
  }

  for (auto &vp : gz.vertex_properties()) {
    vp.second.erase("width");
    vp.second.erase("height");
  }

  /*
  for (V v : g.vertices()) {
    for (auto e : ads.sketches[v]) {
      if (e.d > 1) gz.set_edge_property(v, e.v, "label", to_string(e.d));
    }
  }
  */

  gz.generate_png("02_ads_and_srs.png");

  return 0;
}
