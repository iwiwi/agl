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

    W radius = agl::random(3) + 1;
    const int k = 128;
    vector<V> rank(g.num_vertices());
    vector<V> inv(g.num_vertices());
    for (V i = 0; i < g.num_vertices(); ++i) {
      inv[i] = i;
    }
    random_shuffle(inv.begin(), inv.end());
    for (int i = 0; i < g.num_vertices(); ++i) {
      rank[inv[i]] = i;
    }

    coverage_manager cm(g, radius, 1.0);
    for (int cover_trial = 0; cover_trial < 10; ++cover_trial) {
      vector<bool> covered(g.num_vertices());
      for (int i = 0; i < g.num_vertices(); ++i) covered[i] = cm.v_covered(i);

      vector<vector<V>> naive_x =
          naive_build_sketch(g, radius, k, rank, inv, covered);
      vector<vector<V>> x = build_sketch(g, radius, k, rank, inv, cm);
      for (V v = 0; v < g.num_vertices(); v++) {
        ASSERT_EQ(naive_x[v], x[v]) << v;
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

    coverage_manager cm(g, radius, 1.0);
    vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, cm);

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
      coverage_manager cm(g, radius, 1.0);
      select_greedily(g, X, centers1, centered, k, cm);
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

    coverage_manager cm(g, radius, 1.0);
    vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, cm);
    // pretty_print(g);
    cerr << "radius: " << radius << endl;
    vector<V> centers1;
    {
      double timer = -get_current_time_sec();
      vector<bool> centered(g.num_vertices(), false);
      coverage_manager cm(g, radius, 1.0);
      select_greedily(g, X, centers1, centered, k, cm);
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
    cerr << "size: " << centers1.size() << endl;
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

  coverage_manager cm(g, radius, 1.0);
  vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, cm);
  pretty_print(g);
  cerr << "radius: " << radius << endl;
  vector<V> centers1;
  {
    double timer = -get_current_time_sec();
    vector<bool> centered(g.num_vertices(), false);
    coverage_manager cm(g, radius, 1.0);
    select_greedily(g, X, centers1, centered, k, cm);
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
  cerr << "size: " << centers1.size() << endl;

  ASSERT_EQ(centers1, centers2);
}

TEST(box_cover, coverage_management) {
  for (int trial = 0; trial < 100; ++trial) {
    const int k = 1024;
    W radius = agl::random(2) + 1;
    V M = 3;
    V N = M + agl::random(1000) + 1000;
    while (N <= k) N = M + agl::random(1000) + 1000;
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

    coverage_manager cmtmp(g, radius, 1.0);
    vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, cmtmp);
    pretty_print(g);
    cerr << "radius: " << radius << endl;
    vector<V> centers;
    vector<bool> centered(g.num_vertices(), false);
    coverage_manager cm(g, radius, 1.0);
    select_greedily(g, X, centers, centered, k, cm);
    double tester = coverage(g, centers, radius);
    double hey = cm.get_current_coverage();
    ASSERT_EQ(tester, hey);
  }
}

TEST(box_cover, coverage_break) {
  for (int trial = 0; trial < 100; ++trial) {
    const int k = 1024;
    W radius = agl::random(2) + 1;
    V M = 3;
    V N = M + agl::random(1000) + 1000;
    while (N <= k) N = M + agl::random(1000) + 1000;
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

    coverage_manager cmtmp(g, radius, 1.0);
    vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, cmtmp);
    pretty_print(g);
    vector<bool> centered(g.num_vertices(), false);

    double goal = (double)(agl::random(5) + 95) / 100;
    cerr << "goal: " << goal << endl;
    vector<V> centers = box_cover_sketch(g, radius, k, 100, goal);

    coverage_manager cm(g, radius, goal);
    for (V c : centers) cm.add(g, c);
    ASSERT_TRUE(cm.get_current_coverage() >= goal);
  }
}

TEST(box_cover, solution_flower) {
  vector<pair<V, V>> uvs = {{1, 2}, {2, 2}};
  for (auto p : uvs) {
    V u = p.first, v = p.second;
    V req = 10000;
    auto es = generate_uv_flower(req, u, v);
    G g(make_undirected(es));
    auto pairs = find_analytical_solution("flower", u, v, g);
    for (auto p : pairs) {
      cerr << p.first << " " << p.second << endl;
    }
  }
}

TEST(box_cover, solution_shm) {
  vector<pair<V, int>> nts = {{5, 2}};
  for (auto p : nts) {
    V n = p.first, t = p.second;
    V req = 50;
    auto es = generate_shm(req, n, t);
    G g(make_undirected(es));
    auto pairs = find_analytical_solution("shm", n, t, g);
    for (auto p : pairs) {
      cerr << p.first << " " << p.second << endl;
    }
  }
}

TEST(box_cover, lazy_greedily) {
  V u = 2, v = 2;
  auto es = generate_uv_flower(10000, u, v);
  G g(make_undirected(es));
  vector<pair<W, V>> pairs = find_analytical_solution("flower", u, v, g);

  vector<V> rank(g.num_vertices());
  vector<V> inv(g.num_vertices());
  for (V i = 0; i < g.num_vertices(); ++i) {
    inv[i] = i;
  }
  random_shuffle(inv.begin(), inv.end());
  for (int i = 0; i < g.num_vertices(); ++i) {
    rank[inv[i]] = i;
  }
  cerr << g.num_vertices() << endl;
  double timer = 0;
  const int k = 1024;
  for (auto p : pairs) {
    W rad = p.first / 2;
    cerr << rad << endl;

    timer = -get_current_time_sec();
    coverage_manager cm(g, rad, 0.98);
    auto X = build_sketch(g, rad, k, rank, inv, cm);
    timer += get_current_time_sec();
    cerr << timer << endl;
    timer = -get_current_time_sec();
    double cnt = 0.0;
    for (auto x : X)
      if (x.size() == k) cnt += 1.0;
    cnt /= X.size();

    vector<V> rank(g.num_vertices());
    vector<V> inv(g.num_vertices());
    for (V i = 0; i < g.num_vertices(); ++i) {
      inv[i] = i;
    }
    random_shuffle(inv.begin(), inv.end());
    for (int i = 0; i < g.num_vertices(); ++i) {
      rank[inv[i]] = i;
    }
    vector<V> centers;
    select_lazy_greedily(g, X, rank, inv, centers, cm);
    timer += get_current_time_sec();
    cerr << cnt << endl;
    cerr << centers.size() << " " << timer << endl;

    timer = -get_current_time_sec();
    auto memb = box_cover_memb(g, rad);
    timer += get_current_time_sec();
    cerr << memb.size() << " " << timer << endl;

    cerr << p.second << endl;
  }
}

TEST(box_cover, coloring) {
  V u = 2, v = 2;
  auto es = generate_uv_flower(1000, u, v);
  G g(make_undirected(es));
  pretty_print(g);

  vector<pair<W, V>> pairs = find_analytical_solution("flower", u, v, g);
  for (auto p : pairs) {
    cerr << p << endl;
  }

  W largest = pairs.rbegin()->first;
  vector<pair<W, size_t>> coloring = box_cover_coloring(g, largest);
  for (auto c : coloring) {
    cerr << c << endl;
  }
}

TEST(box_cover, fast_build_sketch) {
  for (int trial = 0; trial < 5; ++trial) {
    const W radius = 256;
    const int k = 128;
    V u = 2, v = 2;
    auto es = generate_uv_flower(10000, u, v);
    G g(make_undirected(es));
    pretty_print(g);

    vector<V> rank(g.num_vertices());
    vector<V> inv(g.num_vertices());
    for (V i = 0; i < g.num_vertices(); ++i) {
      inv[i] = i;
    }
    random_shuffle(inv.begin(), inv.end());
    for (int i = 0; i < g.num_vertices(); ++i) {
      rank[inv[i]] = i;
    }
    coverage_manager cm(g, radius, 1.0);

    cerr << "k: " << k << endl;
    cerr << "radius: " << radius << endl;
    double timer = -get_current_time_sec();
    vector<vector<V>> bX = build_sketch(g, radius, k, rank, inv, cm);
    timer += get_current_time_sec();
    cerr << "build_sketch():" << timer << " sec" << endl;
    bool use_memb = false;

    timer = -get_current_time_sec();
    vector<vector<V>> fX =
        fast_build_sketch(g, radius, k, rank, cm, use_memb, 1000000000);
    timer += get_current_time_sec();
    cerr << "fast_build_sketch():" << timer << " sec" << endl;

    if (bX != fX) {
      int N = bX.size();
      for (int i = 0; i < N; ++i) {
        if (bX[i].size() != fX[i].size()) {
          cerr << i << endl;
          cerr << bX[i].size() << " " << fX[i].size() << endl;
          V mi = min(bX[i].size(), fX[i].size());
          for (int v = 0; v < mi; ++v) {
            cerr << bX[i][v] << " " << fX[i][v] << endl;
          }
        }
      }
    }
    ASSERT_EQ(bX, fX);
  }
}

TEST(box_cover, covered_check) {
  for (int trial = 0; trial < 10; ++trial) {
    V N = agl::random(1000) + 1000;
    W radius = agl::random(10) + 10;
    int k = 128;
    auto es = generate_uv_flower(N, 2, 2);
    G g(make_undirected(es));
    coverage_manager cm(g, radius, 1.0);
    vector<V> rank(g.num_vertices());
    vector<V> inv(g.num_vertices());
    for (V i = 0; i < g.num_vertices(); ++i) inv[i] = i;
    random_shuffle(inv.begin(), inv.end());
    for (int i = 0; i < g.num_vertices(); ++i) rank[inv[i]] = i;

    cm.add(g, 0);
    vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, cm);
    ASSERT_TRUE(X[0].empty());
  }
}