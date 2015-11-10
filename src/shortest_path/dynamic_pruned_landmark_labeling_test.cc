// typed test を使って複数のテンプレートパラメータでテストして下さい
// 参考：https://github.com/iwiwi/pruned-landmark-labeling/blob/master/src/pruned_landmark_labeling_test.cc

// 更新のワークロード生成は |tutorial/07_dynamic_update_benchmark_main.cc|
// を参考に

#include <gtest/gtest.h>
#include "dynamic_pruned_landmark_labeling.h"

using namespace std;
using namespace agl;
using testing::Types;

vector<vector<W>> warshall_floyd(const G &g, W INF) {
  V n = g.num_vertices();
  const auto es = g.edge_list();
  vector<vector<W>> dist(n, vector<W>(n, INF));
  for (int i = 0; i < (int)es.size(); ++i) {
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

TEST(dpll_static, query_ba) {
  for (int trial = 0; trial < 1000; ++trial) {
    V M = 3 + agl::random(10);
    V N = agl::random(10) + M;
    unweighted_edge_list es = generate_ba(N, M);
    dynamic_pruned_landmark_labeling<0> dpll;
    G g(es);
    dpll.construct(g);

    vector<vector<W>> dist = warshall_floyd(g, 100);
    V n = g.num_vertices();
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        W dd = dpll.query_distance(g, i, j);
        ASSERT_EQ(dist[i][j], dd) << i << " " << j;
      }
    }
  }
}

typedef Types<dynamic_pruned_landmark_labeling<0>> dpll_types;

template <typename T>
class dpll_test : public testing::Test {};
TYPED_TEST_CASE(dpll_test, dpll_types);

template <typename TypeParam>
void Test(const unweighted_edge_list &es, bool &check) {
  // CONSTRUCTION
  G g(es);
  TypeParam dpll;
  dpll.construct(g);
  vector<vector<W>> dist = warshall_floyd(g, 100);
  V n = g.num_vertices();
  for (V i = 0; i < n; ++i) {
    for (V j = 0; j < n; ++j) {
      W qd = dpll.query_distance(g, i, j);
      W wd = dist[i][j];
      if (qd != wd) {
        auto ess = g.edge_list();
        for (auto e : ess) {
          cerr << "{" << e.first << "," << e.second << "}," << endl;
        }

        for (int ti = 0; ti < n; ++ti) {
          for (int tj = 0; tj < n; ++tj) {
            if (dist[ti][tj] != dpll.query_distance(g, ti, tj)) {
              cerr << ti << "->" << tj << endl;
              cerr << "dpll:	" << dpll.query_distance(g, ti, tj) << endl;
              cerr << "WF:	" << dist[ti][tj] << endl;
            }
          }
        }
        check = false;
        return;
      }
    }
  }

  // ONLINE UPDATE
  for (int q = 0; q < 20; ++q) {
    V q_from = -1, q_to = -1;
    int cnt = 0;
    while (q_from == q_to || dpll.query_distance(g, q_from, q_to) == 1) {
      q_from = agl::random(n);
      q_to = agl::random(n);

      cnt++;
      if (cnt == 100) {
        return;
      }
    }

    W before = dpll.query_distance(g, q_from, q_to);
    dpll.add_edge(g, q_from, q_to);
    update_warshall_floyd(dist, q_from, q_to);

    for (V i = 0; i < n; ++i) {
      for (V j = 0; j < n; ++j) {
        W qd = dpll.query_distance(g, i, j);
        W wd = dist[i][j];
        if (wd != qd) {
          auto ess = g.edge_list();
          for (auto e : ess) {
            cerr << "{" << e.first << "," << e.second << "}," << endl;
          }

          cerr << i << "->" << j << endl;
          cerr << "Add:	" << q_from << "->" << q_to << endl;
          cerr << "before:	" << before << endl;
          cerr << "dpll:	" << qd << endl;
          cerr << "WF:	" << wd << endl;
          check = false;
          return;
        }
      }
    }
  }
  check = true;
}

TYPED_TEST(dpll_test, three_vertices) {
  unweighted_edge_list es;
  es.emplace_back(0, 1);
  es.emplace_back(1, 2);
  G g(es);
  TypeParam dpll;
  dpll.construct(g);
  ASSERT_EQ(dpll.query_distance(g, 0, 2), 2);
  dpll.add_edge(g, 2, 0);
  ASSERT_EQ(dpll.query_distance(g, 2, 0), 1);
}

TYPED_TEST(dpll_test, path) {
  for (int trial = 0; trial < 10; ++trial) {
    V n = agl::random(50) + 5;
    unweighted_edge_list es;
    for (V i = 0; i + 1 < n; ++i) {
      es.emplace_back(i, i + 1);
    }
    bool check = true;
    Test<TypeParam>(es, check);
    ASSERT_TRUE(check);
  }
}

// Circle (0 -- 1 -- 2 --...-- |L-2| -- |L-1| -- 0)
TYPED_TEST(dpll_test, cycle) {
  for (int trial = 0; trial < 10; ++trial) {
    V n = agl::random(50) + 5;
    unweighted_edge_list es = generate_cycle(n);
    bool check = true;
    Test<TypeParam>(es, check);
    ASSERT_TRUE(check);
  }
}

TYPED_TEST(dpll_test, grid) {
  const int N = 10;
  unweighted_edge_list es = generate_grid(N, N);
  G g(es);
  bool check = true;
  Test<TypeParam>(es, check);
  ASSERT_TRUE(check);
}

TYPED_TEST(dpll_test, small_ba) {
  V N = 14;
  V M = 3;
  for (int trial = 0; trial < 1000; ++trial) {
    unweighted_edge_list es = generate_ba(N, M);
    G g(es);
    bool check = true;
    Test<TypeParam>(es, check);
    ASSERT_TRUE(check);
  }
}

TYPED_TEST(dpll_test, random_ba) {
  for (int trial = 0; trial < 1000; ++trial) {
    V M = 3 + agl::random(10);
    V N = agl::random(10) + M;
    unweighted_edge_list es = generate_ba(N, M);
    G g(es);
    bool check = true;
    Test<TypeParam>(es, check);
    ASSERT_TRUE(check);
  }
}

TYPED_TEST(dpll_test, random_hk) {
  for (int trial = 0; trial < 1000; ++trial) {
    V M = 3 + agl::random(10);
    V N = agl::random(10) + M + 1;
    double P = agl::random(1000) / 1000.0;
    unweighted_edge_list es = generate_hk(N, M, P);
    G g(es);
    bool check = true;
    Test<TypeParam>(es, check);
    ASSERT_TRUE(check);
  }
}