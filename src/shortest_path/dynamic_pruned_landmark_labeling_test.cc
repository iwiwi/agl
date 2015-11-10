// typed test を使って複数のテンプレートパラメータでテストして下さい
// 参考：https://github.com/iwiwi/pruned-landmark-labeling/blob/master/src/pruned_landmark_labeling_test.cc

// 更新のワークロード生成は |tutorial/07_dynamic_update_benchmark_main.cc|
// を参考に

#include <gtest/gtest.h>
#include "dynamic_pruned_landmark_labeling.h"

using namespace std;
using namespace agl;
using testing::Types;

typedef Types<dynamic_pruned_landmark_labeling<0>> dpll_types;

template <typename T>
class dpll_test : public testing::Test {};
TYPED_TEST_CASE(dpll_test, dpll_types);

vector<vector<W>> warshall_floyd(const G &g, W INF) {
  V n = g.num_vertices();
  const auto es = g.edge_list();
  vector<vector<W>> dist(n, vector<W>(n, INF));
  for (int i = 0; i < es.size(); ++i) {
    V v_from = es[i].first;
    V v_to = es[i].second;
    dist[v_from][v_to] = 1;
  }
  for (int i = 0; i < n; ++i) {
    dist[i][i] = 0;
  }
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
      }
    }
  }
  return dist;
}

template <typename TypeParam>
void Test(const G &g) {
  TypeParam dpll;
  dpll.construct(g);
  vector<vector<W>> dist = warshall_floyd(g, dpll.INF8);
  V n = g.num_vertices();
  for (V i = 0; i < n; ++i) {
    for (V j = 0; j < n; ++j) {
      W qd = dpll.query_distance(g, i, j);
      W wd = dist[i][j];
      ASSERT_EQ(wd, qd) << i << " " << j;
    }
  }
}

TYPED_TEST(dpll_test, three_vertices) {
  unweighted_edge_list es;
  es.emplace_back(0, 1);
  es.emplace_back(1, 2);
  G g(es);
  Test<TypeParam>(g);
}

TYPED_TEST(dpll_test, path) {
  const V L = 30;
  unweighted_edge_list es;
  for (V i = 0; i + 1 < L; ++i) {
    es.emplace_back(i, i + 1);
  }
  G g(es);
  Test<TypeParam>(g);
}

// Circle (0 -- 1 -- 2 --...-- |L-2| -- |L-1| -- 0)
TYPED_TEST(dpll_test, cycle) {
  const V L = 30;
  unweighted_edge_list es = generate_cycle(L);
  Test<TypeParam>(G(es));
}

// Almost empty (for testing disconnected pairs)
TYPED_TEST(dpll_test, almost_empty) {
  unweighted_edge_list es;
  es.emplace_back(0, 3);
  es.emplace_back(3, 6);
  es.emplace_back(6, 0);
  Test<TypeParam>(G(es));
}

// // Erdos-Renyi random graph
// TYPED_TEST(dpll_test, famous_model) {
//   const int N = 50;
//   unweighted_edge_list es;
//   for (int v = 0; v < N; ++v) {
//     for (int w = v + 1; w < N; ++w) {
//       double x = rand() / double(RAND_MAX);
//       if (x < 0.5) {
//         es.push_back(make_pair(v, w));
//       }
//     }
//   }
//   Test<TypeParam>(N, es);
// }