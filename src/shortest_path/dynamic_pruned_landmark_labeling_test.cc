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
void cerrIndex(const G& g, TypeParam& dpll) {
  V num_v = g.num_vertices();
  for (int direction = 0; direction < 2; ++direction) {
    cerr << direction << endl;
    for (V v = 0; v < num_v; ++v) {
      cerr << v << " ";
      for (const auto& p : dpll.idx[direction][v].spm) {
        cerr << "[" << p.first << "," << p.second << "] ";
      }
      cerr << endl;
    }
  }
}

template <typename TypeParam>
bool BFSCheck(const G& g, TypeParam& dpll) {
  // pretty_print(g);
  double t = 0;
  const W INF = 100;
  V num_v = g.num_vertices();
  int cnt = 0;
  for (V from = 0; from < num_v; ++from) {
    if (cnt < from * 100 / num_v) {
      cnt++;
      cerr << "#";
    }
    vector<W> dist(num_v, INF);
    queue<V> que;
    que.push(from);
    dist[from] = 0;
    while (!que.empty()) {
      V v = que.front();
      que.pop();
      for (const auto& u : g.neighbors(v)) {
        if (dist[u] < INF) continue;
        dist[u] = dist[v] + 1;
        que.push(u);
      }
    }
    t = -get_current_time_sec();
    for (int j = 0; j < num_v; ++j) {
      if (dist[j] == INF) continue;
      EXPECT_EQ(dist[j], dpll.query_distance(g, from, j)) << from << "->" << j;
      if (dist[j] != dpll.query_distance(g, from, j)) {
        cerrIndex(g, dpll);
        return false;
      }
    }
    t += get_current_time_sec();
  }
  cerr << endl;
  t /= num_v * num_v;
  cerr << t * 1000 * 1000 << " micro sec/query" << endl;
  return true;
}

template <typename TypeParam>
bool Test(const G& g) {
  TypeParam dpll;
  double t = -get_current_time_sec();
  dpll.construct(g);
  t += get_current_time_sec();
  cerr << t << " sec/construct" << endl;

  return BFSCheck(g, dpll);
}

template <typename TypeParam>
bool DynamicTest(const G& g) {
  // pretty_print(g);
  TypeParam dpll;
  double t = -get_current_time_sec();
  dpll.construct(g);
  t += get_current_time_sec();
  cerr << t << " sec/construct" << endl;

  auto es = g.edge_list();
  V trial = min(100, g.num_vertices() / 2);
  t = -get_current_time_sec();
  vector<pair<V, V>> addition;
  for (int add_trial = 0; add_trial < trial; ++add_trial) {
    V v_from = agl::random(g.num_vertices());
    V v_to = agl::random(g.num_vertices());
    dpll.add_edge(g, v_from, v_to);
    es.emplace_back(v_from, v_to);
    addition.emplace_back(v_from, v_to);
  }
  t += get_current_time_sec();
  cerr << (t / trial) << " sec/add" << endl;
  G g_mod(es);

  if (BFSCheck(g_mod, dpll)) {
    return true;
  } else {
    set<pair<V, V>> s;
    for (const auto& p : g.edge_list()) {
      cerr << "{" << p.first << "," << p.second << "},";
      s.insert(p);
    }
    cerr << endl;
    cerr << "---" << endl;
    for (const auto& p : addition) {
      cerr << "{" << p.first << "," << p.second << "},";
    }
    cerr << endl;
    pretty_print(g);
    pretty_print(g_mod);
    return false;
  }
}

TYPED_TEST(dpll_test, small_grid) {
  auto es = generate_grid(3, 3);
  G g(es);
  ASSERT_TRUE(Test<TypeParam>(g));
}

TYPED_TEST(dpll_test, small_ba) {
  for (int trial = 0; trial < 100; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(50) + m;
    auto es = generate_ba(n, m);
    G g(es);
    ASSERT_TRUE(Test<TypeParam>(g));
  }
}

TYPED_TEST(dpll_test, medium_ba) {
  for (int trial = 0; trial < 10; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(1000) + m;
    auto es = generate_ba(n, m);
    G g(es);
    ASSERT_TRUE(Test<TypeParam>(g));
  }
}

TYPED_TEST(dpll_test, large_ba) {
  for (int trial = 0; trial < 5; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(10000) + m;
    auto es = generate_ba(n, m);
    G g(es);
    pretty_print(g);
    ASSERT_TRUE(Test<TypeParam>(g));
  }
}

TYPED_TEST(dpll_test, dynamic_small_grid) {
  auto es = generate_grid(3, 3);
  G g(es);
  TypeParam dpll;
  dpll.construct(g);
  ASSERT_TRUE(DynamicTest<TypeParam>(g));
}

TYPED_TEST(dpll_test, dynamic_tedukuri) {
  unweighted_edge_list es = {
      {0, 1},  {0, 2},  {0, 6}, {1, 2}, {1, 3},  {1, 4},  {1, 7},  {2, 3},
      {2, 4},  {2, 5},  {2, 7}, {2, 8}, {2, 11}, {3, 5},  {3, 9},  {3, 10},
      {3, 11}, {3, 12}, {4, 6}, {6, 9}, {7, 8},  {8, 12}, {9, 10},
  };
  G g1(es);
  TypeParam dpll;
  dpll.construct(g1);

  unweighted_edge_list addition = {
      {5, 10}, {0, 10}, {6, 5}, {10, 8}, {7, 9}, {0, 2},
  };
  for (const auto& p : addition) {
    es.emplace_back(p.first, p.second);
    dpll.add_edge(g1, p.first, p.second);
  }
  G g2(es);
  ASSERT_TRUE(BFSCheck(g2, dpll));
}

TYPED_TEST(dpll_test, dynamic_small_ba) {
  for (int trial = 0; trial < 100; ++trial) {
    V m = 2;
    V n = agl::random(20) + 1 + m;
    auto es = generate_ba(n, m);
    G g(es);
    ASSERT_TRUE(DynamicTest<TypeParam>(g));
  }
}

TYPED_TEST(dpll_test, dynamic_medium_ba) {
  for (int trial = 0; trial < 10; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(1000) + m;
    auto es = generate_ba(n, m);
    G g(es);
    ASSERT_TRUE(DynamicTest<TypeParam>(g));
  }
}

TYPED_TEST(dpll_test, dynamic_large_ba) {
  for (int trial = 0; trial < 5; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(10000) + m;
    auto es = generate_ba(n, m);
    G g(es);
    pretty_print(g);
    ASSERT_TRUE(DynamicTest<TypeParam>(g));
  }
}