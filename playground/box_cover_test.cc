#include "box_cover.h"
#include "gtest/gtest.h"
#include <sys/time.h>

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

TEST(box_cover, greedy_small) {
  for (int trial = 0; trial < 1000; ++trial) {
    const W radius = agl::random(3) + 1;
    V M = 3;
    V N = M + agl::random(50);
    auto es = generate_ba(N, M);
    G g(make_undirected(es));
    const int k = agl::random(20) + 5;
    if (g.num_vertices() <= k) {
      trial--;
      continue;
    }
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

    vector<V> centers2;
    {
      cerr << "naive" << endl;
      vector<bool> centered(g.num_vertices(), false);
      naive_select_greedily(g, X, centers2, centered, k);
      cerr << centers2 << endl;
    }

    vector<V> centers1;
    {
      cerr << "greedy" << endl;
      vector<bool> centered(g.num_vertices(), false);
      select_greedily(g, X, centers1, centered, k);
    }
    ASSERT_EQ(centers1, centers2);
    cerr << "Stage: " << (trial + 1) << " CLEARED!!!!!!!" << endl;
  }
}

TEST(box_cover, greedy_big) {
  for (int trial = 0; trial < 30; ++trial) {
    const W radius = agl::random(4) + 1;
    V M = 3;
    V N = M + agl::random(2000) + 2000;
    auto es = generate_ba(N, M);
    G g(make_undirected(es));
    const int k = 1024;
    if (g.num_vertices() < k) {
      trial--;
      continue;
    }
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
    pretty_print(g);
    cerr << "radius: " << radius << endl;
    vector<V> centers1;
    {
      double timer = -get_current_time_sec();
      vector<bool> centered(g.num_vertices(), false);
      select_greedily(g, X, centers1, centered, k);
      timer += get_current_time_sec();
      cerr << "greedy: " << timer << " sec" << endl;
    }
    vector<V> centers2;
    {
      double timer = -get_current_time_sec();
      vector<bool> centered(g.num_vertices(), false);
      naive_select_greedily(g, X, centers2, centered, k);
      timer += get_current_time_sec();
      cerr << "naive: " << timer << " sec" << endl;
    }

    ASSERT_EQ(centers1, centers2);
    cerr << "Stage: " << (trial + 1) << " CLEARED!!!!!!!" << endl;
  }
}

TEST(box_cover, greedy_huge) {
  const int k = 1024;
  W radius = agl::random(4) + 1;
  V M = 3;
  V N = M + agl::random(10000) + 10000;
  while (N <= k) N = M + agl::random(10000) + 10000;
  auto es = generate_ba(N, M);
  G g(make_undirected(es));
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
  pretty_print(g);
  cerr << "radius: " << radius << endl;
  vector<V> centers1;
  {
    double timer = -get_current_time_sec();
    vector<bool> centered(g.num_vertices(), false);
    select_greedily(g, X, centers1, centered, k);
    timer += get_current_time_sec();
    cerr << "greedy: " << timer << " sec" << endl;
  }
  vector<V> centers2;
  {
    double timer = -get_current_time_sec();
    vector<bool> centered(g.num_vertices(), false);
    naive_select_greedily(g, X, centers2, centered, k);
    timer += get_current_time_sec();
    cerr << "naive: " << timer << " sec" << endl;
  }

  ASSERT_EQ(centers1, centers2);
}