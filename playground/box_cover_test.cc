#include "box_cover.h"
#include "gtest/gtest.h"
using namespace agl;
using namespace std;

TEST(box_cover, memb) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = 3;
    V N = M + agl::random(1000);
    auto es = generate_ba(N, M);
    G g(make_undirected(es));
    pretty_print(g);

    W radius = 3;
    vector<V> memb = box_cover_memb(g, radius);

    ASSERT_EQ(coverage(g, memb, radius), 1.0);
  }
}

TEST(box_cover, burning) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = 3;
    V N = M + agl::random(200);
    auto es = generate_ba(N, M);
    G g(make_undirected(es));
    pretty_print(g);

    W radius = 1;
    vector<V> burning = box_cover_burning(g, radius);
    vector<W> central_distances(g.num_vertices(), g.num_vertices());

    ASSERT_EQ(coverage(g, burning, radius), 1.0);
  }
}

TEST(box_cover, build_sketch_check) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = 3;
    V N = M + agl::random(1000);
    auto es = generate_ba(N, M);
    G g(make_undirected(es));
    pretty_print(g);

    W radius = 1;
    const int k = 200;
    vector<V> rank(g.num_vertices());
    vector<V> inv(g.num_vertices());
    for (V i = 0; i < g.num_vertices(); ++i) {
      inv[i] = i;
    }
    random_shuffle(inv.begin(), inv.end());
    for (int i = 0; i < g.num_vertices(); ++i) {
      rank[inv[i]] = i;
    }

    vector<bool> covered(g.num_vertices(), false);
    for (int cover_trial = 0; cover_trial < 10; ++cover_trial) {
      vector<vector<V>> naive_x =
          naive_build_sketch(g, radius, k, rank, inv, covered);
      vector<vector<V>> x = build_sketch(g, radius, k, rank, inv, covered);
      for (V v = 0; v < g.num_vertices(); v++) {
        ASSERT_EQ(naive_x[v], x[v]) << v;
      }

      for (int i = 0; i < 100; ++i) {
        int target = agl::random(g.num_vertices());
        covered[target] = true;
      }
    }
  }
}

TEST(box_cover, small_select_greedily) {
  const W radius = 1;
  auto es = generate_grid(4, 5);
  G g(make_undirected(es));
  const int k = 3;
  vector<V> rank(g.num_vertices());
  vector<V> inv(g.num_vertices());
  for (V i = 0; i < g.num_vertices(); ++i) {
    inv[i] = i;
  }
  random_shuffle(inv.begin(), inv.end());
  for (int i = 0; i < g.num_vertices(); ++i) {
    rank[inv[i]] = i;
  }
  vector<bool> covered(g.num_vertices(), false);
  vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, covered);
  for (int i = 0; i < g.num_vertices(); ++i) {
    cerr << i << " " << X[i] << endl;
  }

  vector<V> centers1;
  {
    vector<bool> centered(g.num_vertices(), false);
    select_greedily(g, X, inv, centers1, centered, k);
  }
  vector<V> centers2;
  {
    vector<bool> centered(g.num_vertices(), false);
    naive_select_greedily(g, X, centers2, centered, k);
  }
  cerr << centers1 << endl;
  cerr << centers2 << endl;
  ASSERT_EQ(centers1, centers2);
}

TEST(box_cover, large_select_greedily) {
  const W radius = 3;
  auto es = generate_grid(30, 40);
  G g(make_undirected(es));
  const int k = 50;
  vector<V> rank(g.num_vertices());
  vector<V> inv(g.num_vertices());
  for (V i = 0; i < g.num_vertices(); ++i) {
    inv[i] = i;
  }
  random_shuffle(inv.begin(), inv.end());
  for (int i = 0; i < g.num_vertices(); ++i) {
    rank[inv[i]] = i;
  }
  vector<bool> covered(g.num_vertices(), false);
  vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, covered);

  vector<V> centers1;
  {
    vector<bool> centered(g.num_vertices(), false);
    select_greedily(g, X, inv, centers1, centered, k);
  }
  vector<V> centers2;
  {
    vector<bool> centered(g.num_vertices(), false);
    naive_select_greedily(g, X, centers2, centered, k);
  }
  cerr << centers1 << endl;
  cerr << centers2 << endl;
  ASSERT_EQ(centers1, centers2);
}
