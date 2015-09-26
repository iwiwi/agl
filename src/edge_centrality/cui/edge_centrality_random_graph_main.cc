#include "easy_cui.h"
#include "edge_centrality/edge_centrality.h"

DEFINE_int32(n, 10, "n");
DEFINE_int64(m, 30, "m");

namespace {
unweighted_edge_list gen_connected_random_graph(V num_vs, size_t num_es) {
  CHECK(num_es + 1 >= (size_t)num_vs);
  CHECK(num_es <= num_vs * (size_t)(num_vs - 1));

  unordered_set<pair<V, V>> es;

  // Generate a spanning tree
  unweighted_edge_list spanning_es = generate_random_spanning_tree(num_vs);
  for (auto e : spanning_es) {
    es.emplace(min(e.first, e.second), max(e.first, e.second));
  }

  // Add random edges
  std::uniform_int_distribution<V> rng(0, num_vs - 1);
  while (es.size() < num_es) {
    V u = rng(agl::random), v = rng(agl::random);
    if (u > v) swap(u, v);
    if (u != v) es.insert(make_pair(u, v));
  }
  return vector<pair<V, V>>(es.begin(), es.end());
}
}  // namespace

int main(int argc, char **argv) {
  JLOG_INIT(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);

  G g(make_undirected(gen_connected_random_graph(FLAGS_n, FLAGS_m)));
  // CHECK(is_connected(g));

  auto ebc = edge_betweenness_centrality(g);
  ebc = merge_edge_centrality_map_entries_for_undirected_graph(ebc);

  for (auto i : ebc) {
    cout << i.first.first << "\t" << i.first.second << "\t" << i.second << endl;
  }

  return 0;
}
