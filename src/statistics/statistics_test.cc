#include "statistics.h"
#include "gtest/gtest.h"
#include "prettyprint.h"
#include "visualize/graphviz.h"
using namespace agl;
using namespace std;

namespace {
pair<size_t, size_t> num_triangles_and_wedges_naive(const G &g) {
  size_t num_triangles = 0;
  size_t num_wedges = 0;

  for (V x : g.vertices()) {
    for (V y : g.vertices()) {
      for (V z : g.vertices()) {
        if (x < y && y < z) {
          int c = 0;
          c += (is_adjacent(g, x, y) || is_adjacent(g, y, x));
          c += (is_adjacent(g, y, z) || is_adjacent(g, z, y));
          c += (is_adjacent(g, z, x) || is_adjacent(g, x, z));
          if (c == 3) ++num_triangles;
          if (c == 2) ++num_wedges;
        }
      }
    }
  }
  return make_pair(num_triangles, num_wedges);
}

vector<pair<size_t, size_t>> num_local_triangles_and_wedges_naive(const G &g) {
  vector<pair<size_t, size_t>> res(g.num_vertices());

  for (V x : g.vertices()) {
    for (V y : g.vertices()) {
      for (V z : g.vertices()) {
        if (x < y && y < z) {
          int c = 0;
          c += (is_adjacent(g, x, y) || is_adjacent(g, y, x));
          c += (is_adjacent(g, y, z) || is_adjacent(g, z, y));
          c += (is_adjacent(g, z, x) || is_adjacent(g, x, z));
          if (c == 3) {
            ++res[x].first;
            ++res[y].first;
            ++res[z].first;
          }
          if (c == 2) {
            if (!(is_adjacent(g, x, y) || is_adjacent(g, y, x))) ++res[z].second;
            if (!(is_adjacent(g, y, z) || is_adjacent(g, z, y))) ++res[x].second;
            if (!(is_adjacent(g, z, x) || is_adjacent(g, x, z))) ++res[y].second;
          }
        }
      }
    }
  }
  return res;
}
}  // namespace

TEST(num_triangles_and_wedges, random) {
  for (int trial = 0; trial < 10; ++trial) {
    G g(generate_erdos_renyi(100, 5));
    ASSERT_EQ(num_triangles_and_wedges_naive(g), num_triangles_and_wedges(g));
  }
}

TEST(num_local_triangles_and_wedges, random) {
  for (int trial = 0; trial < 10; ++trial) {
    G g(generate_erdos_renyi(100, 5));
    ASSERT_EQ(num_local_triangles_and_wedges_naive(g), num_local_triangles_and_wedges(g));
  }
}
