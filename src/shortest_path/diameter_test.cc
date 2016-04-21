#include "diameter.h"
#include "dijkstra.h"
#include "gtest/gtest.h"
using namespace std;
using namespace agl;

int dijkstra_diameter(const G& g) {
  int diameter = 0;
  auto dm = all_pairs_distance(g);
  for (const auto& v : g.vertices()) {
    for (V u : irange<V>(g.num_vertices())) {
      if (dm[v][u] != numeric_limits<int>::max()) {
        diameter = max(diameter, dm[v][u]);
      }
    }
  }
  return diameter;
}

TEST(diameter, grid) {
  auto es = generate_grid(3, 3);
  G g(es);
  ASSERT_EQ(diameter(g), 4);
}

TEST(diameter, ba) {
  for (int trial = 0; trial < 100; ++trial) {
    V m = agl::random(10) + 2;
    V n = agl::random(50) + m + 1;
    auto es = generate_ba(n, m);
    G g(es);
    ASSERT_EQ(diameter(g), dijkstra_diameter(g));
  }
}
