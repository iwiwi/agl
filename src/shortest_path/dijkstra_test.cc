#include "dijkstra.h"
#include "gtest/gtest.h"
using namespace std;
using namespace agl;

TEST(dijkstra, unweighted) {
  G g(G::edge_list_type({
          {0, {1}},
          {1, {2}},
          {2, {3}},
          {3, {1}},
  }));

  auto ds = single_source_distance(g, 0);
  for (auto d : ds) {
    cout << d << " ";
  }
  cout << endl;
}

TEST(dijkstra, weighted) {
  weighted_graph<> g(weighted_graph<>::edge_list_type{
          {0, {1, 1.5}},
          {1, {2, 20}},
          {2, {3, 5}},
          {3, {1, 0.01}},
  });

  auto ds = single_source_distance(g, 0);
  for (auto d : ds) {
    cout << d << " ";
  }
  cout << endl;
}

namespace {
int fac(int k) {
  return k == 0 ? 1 : k * fac(k - 1);
}

int cmb(int n, int k) {
  return fac(n) / fac(k) / fac(n - k);
}
}

TEST(single_source_distance_with_num_paths, grid) {
  static constexpr int kNumRows = 5, kNumCols = 5;
  auto vid = [=](int i, int j) { return i * kNumCols + j; };

  G g(gen_grid(kNumRows, kNumCols));
  auto r = single_source_distance_with_num_paths(g, 0);

  for (auto i : make_irange(kNumRows)) {
    for (auto j : make_irange(kNumCols)) {
      int v = vid(i, j);
      ASSERT_EQ(r[v].first, i + j);
      ASSERT_EQ(r[v].second, cmb(i + j, i));
    }
  }
}
