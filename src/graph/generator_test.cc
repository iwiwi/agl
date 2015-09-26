#include "generator.h"
#include "connectivity/connectivity.h"
#include "gtest/gtest.h"
using namespace agl;

TEST(gen_radom_spanning_tree, connectivity) {
  for (int trial : make_irange(10)) {
    V num_vs = 1 + agl::random(100000);
    auto es = generate_random_spanning_tree(num_vs);

    // Number of edges
    ASSERT_EQ((V)es.size(), num_vs - 1);

    // Connectivity
    G g(force_undirected(es));
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));
  }
}
