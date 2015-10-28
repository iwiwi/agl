#include "generator.h"
#include "connectivity/connectivity.h"
#include "gtest/gtest.h"
using namespace agl;

TEST(gen_radom_spanning_tree, connectivity) {
  for (int trial = 0; trial < 10; ++trial) {
    V num_vs = 1 + agl::random(100000);
    auto es = generate_random_spanning_tree(num_vs);

    // Number of edges
    ASSERT_EQ((V)es.size(), num_vs - 1);

    // Connectivity
    G g(make_undirected(es));
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));
  }
}

TEST(gen_cycle, some_vertices) {
  auto l0 = generate_cycle(0);
  ASSERT_EQ((V)l0.size(), 0);

  auto l1 = generate_cycle(1);
  ASSERT_EQ((V)l1.size(), 0);

  auto l2 = generate_cycle(2);
  ASSERT_EQ((V)l2.size(), 2);
}

TEST(gen_cycle, random_num_vertices) {
  for (int trial = 0; trial < 10; ++trial) {
    V num_vs = 1 + agl::random(100000);
    auto el = generate_cycle(num_vs);

    // Number of edges
    ASSERT_EQ((V)el.size(), num_vs);

    // Connectivity
    G g(el);
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));
  }
}

TEST(gen_ba, random_num_vertices) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = agl::random(1000);
    V N = M + agl::random(1000);
    auto es = generate_ba(N, M);

    // Number of edges
    ASSERT_EQ((V)es.size(), M * (M - 1) / 2 + (N - M) * M);

    G g(es);
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));
  }
}
