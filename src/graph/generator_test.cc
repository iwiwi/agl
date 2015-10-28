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

TEST(gen_ba, some_vertices) {
  auto es = generate_ba(100, 2);
  G g(es);
  pretty_print(g);
}
