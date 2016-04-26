#include "diameter.h"
#include "gtest/gtest.h"
using namespace std;
using namespace agl;

template<typename GraphT = G>
typename GraphT::W dijkstra_diameter(const GraphT& g) {
  using W = typename GraphT::W;
  W diameter = 0;
  auto dm = all_pairs_distance(g);
  for (const auto& v : g.vertices()) {
    for (V u : irange<V>(g.num_vertices())) {
      if (dm[v][u] != numeric_limits<W>::max()) {
        diameter = max(diameter, dm[v][u]);
      }
    }
  }
  return diameter;
}

TEST(diameter, unweighted_grid) {
  auto es = generate_grid(3, 3);
  G g(es);
  ASSERT_EQ(diameter(g), 4);
}

TEST(diameter, unweighted_ba) {
  for (int trial = 0; trial < 100; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(50) + m + 1;
    auto es = generate_ba(n, m);
    G g(es);
    ASSERT_EQ(diameter(g), dijkstra_diameter(g));
  }
}

TEST(diameter, weighted_grid) {
  using GraphT = weighted_graph<double>;
  auto es = generate_grid(3, 3);
  GraphT g(add_random_weight<GraphT>(es));
  ASSERT_TRUE(is_eq(diameter(g), dijkstra_diameter(g)));
}

TEST(diameter, weighted_ba) {
  using GraphT = weighted_graph<double>;
  for (int trial = 0; trial < 100; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(50) + m + 1;
    auto es = generate_ba(n, m);
    GraphT g(add_random_weight<GraphT>(es));
    ASSERT_TRUE(is_eq(diameter(g), dijkstra_diameter(g)));
  }
}
