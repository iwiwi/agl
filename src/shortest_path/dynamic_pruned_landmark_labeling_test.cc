#include <gtest/gtest.h>
#include "dynamic_pruned_landmark_labeling.h"

using namespace std;
using namespace agl;
using testing::Types;

typedef Types<dynamic_pruned_landmark_labeling<0>> dpll_types;

template <typename T>
class dpll_test : public testing::Test {};
TYPED_TEST_CASE(dpll_test, dpll_types);

template <typename TypeParam>
void Test(const G& g) {
  pretty_print(g);
  TypeParam dpll;
  dpll.construct(g);

  const W INF = 100;
  V num_v = g.num_vertices();
  for (V i = 0; i < num_v; ++i) {
    vector<W> dist(num_v, INF);
    queue<V> que;
    que.push(i);
    dist[i] = 0;
    while (!que.empty()) {
      V v = que.front();
      que.pop();
      for (const auto& u : g.neighbors(v)) {
        if (dist[u] < INF) continue;
        dist[u] = dist[v] + 1;
        que.push(u);
      }
    }
    for (int j = 0; j < num_v; ++j) {
      if (dist[j] == INF) continue;
      cerr << dpll.query_distance(g, i, j) << " " << dist[j] << endl;
      ASSERT_EQ(dist[j], dpll.query_distance(g, i, j)) << i << "->" << j;
    }
  }
}

TYPED_TEST(dpll_test, tekito) {
  auto es = generate_grid(3, 3);
  G g(es);
  Test<TypeParam>(g);
}