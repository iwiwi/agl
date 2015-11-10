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

void update_warshall_floyd(vector<vector<W>> &dist, V v_from, V v_to) {
  V n = dist.size();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      dist[i][j] = min(dist[i][j], dist[i][v_from] + dist[v_to][j] + 1);
    }
  }
}

template <typename TypeParam>
void Test(const unweighted_edge_list &es) {
  G g(es);
  TypeParam dpll;
  dpll.construct(g);
  vector<vector<W>> dist = warshall_floyd(g, dpll.INF8);
  V n = g.num_vertices();
  for (V i = 0; i < n; ++i) {
    for (V j = 0; j < n; ++j) {
      W qd = dpll.query_distance(g, i, j);
      W wd = dist[i][j];
      if (qd != wd) {
        auto ess = g.edge_list();
        for (auto e : ess)
          cerr << "{" << e.first << "," << e.second << "}," << endl;
      }
      ASSERT_EQ(wd, qd) << i << " " << j;
    }
  }

  // for (int q = 0; q < 20; ++q) {
  //   V q_from = -1, q_to = -1;
  //   int cnt = 0;
  //   while (q_from == q_to || dpll.query_distance(g, q_from, q_to) == 1) {
  //     q_from = agl::random(n);
  //     q_to = agl::random(n);

  //     cnt++;
  //     if (cnt == 100) {
  //       return;
  //     }
  //   }

  //   cout << q_from << "->" << q_to << endl;
  //   dpll.add_edge(g, q_from, q_to);
  //   update_warshall_floyd(dist, q_from, q_to);

  //   for (V i = 0; i < n; ++i) {
  //     for (V j = 0; j < n; ++j) {
  //       W qd = dpll.query_distance(g, i, j);
  //       W wd = dist[i][j];
  //       ASSERT_EQ(wd, qd) << i << " " << j;
  //     }
  //   }
  // }
}

TYPED_TEST(dpll_test, three_vertices) {
  unweighted_edge_list es;
  es.emplace_back(0, 1);
  es.emplace_back(1, 2);
  Test<TypeParam>(es);
}

TYPED_TEST(dpll_test, path) {
  for (int trial = 0; trial < 10; ++trial) {
    V n = agl::random(50) + 5;
    unweighted_edge_list es;
    for (V i = 0; i + 1 < n; ++i) {
      es.emplace_back(i, i + 1);
    }
    Test<TypeParam>(es);
  }
}

// Circle (0 -- 1 -- 2 --...-- |L-2| -- |L-1| -- 0)
TYPED_TEST(dpll_test, cycle) {
  for (int trial = 0; trial < 10; ++trial) {
    V n = agl::random(50) + 5;
    unweighted_edge_list es = generate_cycle(n);
    Test<TypeParam>(es);
  }
}

TYPED_TEST(dpll_test, grid) {
  const int N = 10;
  unweighted_edge_list es = generate_grid(N, N);
  G g(es);
  pretty_print(g);
  Test<TypeParam>(es);
}

TYPED_TEST(dpll_test, small_ba) {
  V N = 5;
  V M = 3;
  unweighted_edge_list es = generate_ba(N, M);
  G g(es);
  Test<TypeParam>(es);
}

TYPED_TEST(dpll_test, kanashimi) {
  unweighted_edge_list es = {{0, 1},
                             {0, 2},
                             {0, 3},
                             {0, 4},
                             {0, 8},
                             {0, 9},
                             {0, 11},
                             {1, 2},
                             {1, 3},
                             {1, 4},
                             {1, 6},
                             {1, 11},
                             {1, 13},
                             {2, 3},
                             {2, 4},
                             {2, 5},
                             {2, 8},
                             {2, 12},
                             {2, 13},
                             {3, 5},
                             {3, 8},
                             {3, 10},
                             {3, 12},
                             {3, 14},
                             {4, 5},
                             {4, 6},
                             {4, 7},
                             {4, 10},
                             {4, 12},
                             {5, 6},
                             {5, 7},
                             {5, 9},
                             {5, 14},
                             {6, 7},
                             {6, 11},
                             {7, 9},
                             {7, 10},
                             {7, 14},
                             {8, 13}};
  G g(es);
  pretty_print(g);
  Test<TypeParam>(es);
}
