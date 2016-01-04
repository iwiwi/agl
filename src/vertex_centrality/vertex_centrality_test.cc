#include "vertex_centrality.h"
#include "gtest/gtest.h"
using namespace std;
using namespace agl;

TEST(vertex_centrality, classic_closeness_path) {
  constexpr V kNumVs = 10;
  G g(make_undirected(generate_path(kNumVs)));

  auto cc = vertex_centrality_classic_closeness(g);
  assert(cc.size() == (size_t)kNumVs);

  for (V v : g.vertices()) {
    W s = 0;
    for (V u : g.vertices()) s += abs(v - u);
    ASSERT_NEAR(cc[v], 1.0 / s, 1E-9);
  }
}

TEST(vertex_centrality, pagerank_cycle) {
  constexpr V kNumVs = 10;
  G g(generate_cycle(kNumVs));
  auto pr = vertex_centrality_pagerank(g);
  assert(pr.size() == (size_t)kNumVs);

  for (V v : g.vertices()) {
    ASSERT_NEAR(pr[v], 1.0 / kNumVs, 1E-9);
  }
}

TEST(vertex_centrality, pagerank_skewed) {
  constexpr V kNumVs = 30;
  unweighted_edge_list es;

  for (V v : make_irange(kNumVs)) {
    es.emplace_back(0, v);
    for (V u : make_irange(v)) {
      es.emplace_back(v, u);
    }
  }
  G g(es);
  auto pr = vertex_centrality_pagerank(g);
  assert(pr.size() == (size_t)kNumVs);

  for (V v = 0; v + 1 < kNumVs; ++v) {
    ASSERT_GT(pr[v], pr[v + 1]);
  }
}
