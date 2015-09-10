#include "agl.h"
#include "visitor.h"
#include "gtest/gtest.h"
using namespace std;
using namespace agl;

TEST(visitor_by_distance, unweighted_path) {
  unweighted_edge_list es = gen_path(100);
  G g(es);
  visitor_by_distance<G> vis(g);

  for (int n : {0, 1, 10, 100}) {
    int k = 0;
    vis.visit(0, [&](V v, W w) -> bool {
      // TODO: ASSERT cannot be used in lambda? (probably gtest bug)
      CHECK(k == v);  // ASSERT_EQ(k, v);
      CHECK(v == w);  // ASSERT_EQ(v, w);
      ++k;
      return k <= n;
    });
    ASSERT_EQ(k, min(n + 1, 100));
  }
}

TEST(visitor_by_distance, weighted_path) {
  using graph_type = weighted_graph<double>;
  graph_type::edge_list_type es = add_random_weight<graph_type>(gen_path(100));
  graph_type g(es);
  visitor_by_distance<graph_type> vis(g);

  for (int n : {0, 1, 10, 100}) {
    int k = 0;
    vis.visit(0, [&](V v, graph_type::W w) -> bool {
      CHECK(k == v);  // ASSERT_EQ(k, v);
      ++k;
      return k <= n;
    });
    ASSERT_EQ(k, min(n + 1, 100));
  }
}
