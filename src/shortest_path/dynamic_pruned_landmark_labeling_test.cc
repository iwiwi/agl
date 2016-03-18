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
void BFSCheck(const G& g, TypeParam& dpll) {
  pretty_print(g);
  double t = 0;
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
    t = -get_current_time_sec();
    for (int j = 0; j < num_v; ++j) {
      if (dist[j] == INF) continue;
      ASSERT_EQ(dist[j], dpll.query_distance(g, i, j)) << i << "->" << j;
    }
    t += get_current_time_sec();
  }
  t /= num_v * num_v;
  cerr << t << " sec/query" << endl;
}

template <typename TypeParam>
void Test(const G& g) {
  TypeParam dpll;
  double t = -get_current_time_sec();
  dpll.construct(g);
  t += get_current_time_sec();
  cerr << t << " sec/construct" << endl;

  BFSCheck(g, dpll);
}

template <typename TypeParam>
void DynamicTest(const G& g) {
  pretty_print(g);
  TypeParam dpll;
  double t = -get_current_time_sec();
  dpll.construct(g);
  t += get_current_time_sec();
  cerr << t << " sec/construct" << endl;

  auto es = g.edge_list();
  V trial = 100;
  t = -get_current_time_sec();
  for (int add_trial = 0; add_trial < trial; ++add_trial) {
    V v_from = agl::random(g.num_vertices());
    V v_to = agl::random(g.num_vertices());
    dpll.add_edge(g, v_from, v_to);
    es.emplace_back(v_from, v_to);
  }
  t += get_current_time_sec();
  cerr << (t / trial) << " sec/add" << endl;
  G g_mod(es);

  BFSCheck(g_mod, dpll);
}

TYPED_TEST(dpll_test, small_grid) {
  auto es = generate_grid(3, 3);
  G g(es);
  Test<TypeParam>(g);
}

TYPED_TEST(dpll_test, small_ba) {
  for (int trial = 0; trial < 100; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(50) + m;
    auto es = generate_ba(n, m);
    G g(es);
    Test<TypeParam>(g);
  }
}

TYPED_TEST(dpll_test, medium_ba) {
  for (int trial = 0; trial < 100; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(1000) + m;
    auto es = generate_ba(n, m);
    G g(es);
    Test<TypeParam>(g);
  }
}

TYPED_TEST(dpll_test, large_ba) {
  for (int trial = 0; trial < 5; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(10000) + m;
    auto es = generate_ba(n, m);
    G g(es);
    Test<TypeParam>(g);
  }
}

TYPED_TEST(dpll_test, dynamic_small_grid) {
  auto es = generate_grid(3, 3);
  G g(es);
  TypeParam dpll;
  dpll.construct(g);
  DynamicTest<TypeParam>(g);
}

TYPED_TEST(dpll_test, dynamic_small_ba) {
  for (int trial = 0; trial < 100; ++trial) {
    V m = 3;
    V n = agl::random(7) + m;
    auto es = generate_ba(n, m);
    G g(es);
    DynamicTest<TypeParam>(g);
  }
}

TYPED_TEST(dpll_test, dynamic_medium_ba) {
  for (int trial = 0; trial < 100; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(1000) + m;
    auto es = generate_ba(n, m);
    G g(es);
    DynamicTest<TypeParam>(g);
  }
}

TYPED_TEST(dpll_test, dynamic_large_ba) {
  for (int trial = 0; trial < 5; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(10000) + m;
    auto es = generate_ba(n, m);
    G g(es);
    DynamicTest<TypeParam>(g);
  }
}