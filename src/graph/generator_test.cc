#include "generator.h"
#include "connectivity/connectivity.h"
#include "gtest/gtest.h"
#include <vector>
#include <random>

using namespace agl;
using namespace std;

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

TEST(gen_cycle, some_vertices) {
  auto l0 = generate_cycle(0);
  ASSERT_EQ((V)l0.size(), 0);

  auto l1 = generate_cycle(1);
  ASSERT_EQ((V)l1.size(), 0);

  auto l2 = generate_cycle(2);
  ASSERT_EQ((V)l2.size(), 2);
}

TEST(gen_cycle, random_num_vertices) {
  for (int trial = 0; trial < 10; ++trial) {
    V num_vs = 1 + agl::random(100000);
    auto el = generate_cycle(num_vs);

    // Number of edges
    ASSERT_EQ((V)el.size(), num_vs);

    // Connectivity
    G g(el);
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));
  }
}

TEST(gen_ba, random_num_vertices) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = agl::random(1000) + 3;
    V N = M + agl::random(1000);
    auto es = generate_ba(N, M);

    // Number of edges
    size_t expected_edge_num = M * (M - 1) / 2 + (N - M) * M;
    ASSERT_EQ(es.size(), expected_edge_num);

    G g(es);
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));

    // Check degree
    G ug(make_undirected(es));
    for (V v : ug.vertices()) {
      ASSERT_TRUE(ug.degree(v) >= (size_t)M);
    }
  }
}

TEST(gen_ba, corner_case) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = 2;
    V N = M + agl::random(1000);
    auto es = generate_ba(N, M);

    // Number of edges
    size_t expected_edge_num = M * (M - 1) / 2 + (N - M) * M;
    ASSERT_EQ(es.size(), expected_edge_num);

    G g(es);
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));

    // Check degree
    G ug(make_undirected(es));
    for (V v : ug.vertices()) {
      ASSERT_TRUE(ug.degree(v) >= (size_t)M);
    }
  }
}

TEST(gen_dms, small_case) {
  V N = 5;
  V M = 2;
  V K0 = -1;
  auto es = generate_dms(N, M, K0);

  // Number of edges
  size_t expected_edge_num = M * (M - 1) / 2 + (N - M) * M;
  ASSERT_EQ(es.size(), expected_edge_num);

  G g(es);
  pretty_print(g);
  ASSERT_TRUE(is_connected(g));

  // Check degree
  G ug(make_undirected(es));
  for (V v : ug.vertices()) {
    ASSERT_TRUE(ug.degree(v) >= (size_t)M);
  }
}

TEST(gen_dms, corner_case) {
  V N = 10;
  V M = 0;
  V K0 = -1;
  auto es = generate_dms(N, M, K0);

  // Number of edges
  size_t expected_edge_num = M * (M - 1) / 2 + (N - M) * M;
  ASSERT_EQ(es.size(), expected_edge_num);

  G g(es);
  pretty_print(g);
  ASSERT_TRUE(is_connected(g));

  // Check degree
  G ug(make_undirected(es));
  for (V v : ug.vertices()) {
    ASSERT_TRUE(ug.degree(v) >= (size_t)M);
  }
}

TEST(gen_dms, random_trial) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = agl::random(1000);
    V N = M + agl::random(1000);
    V K0 = M - agl::random(1000);
    auto es = generate_dms(N, M, K0);

    // Number of edges
    size_t expected_edge_num = M * (M - 1) / 2 + (N - M) * M;
    ASSERT_EQ(es.size(), expected_edge_num);

    G g(es);
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));

    // Check degree
    G ug(make_undirected(es));
    for (V v : ug.vertices()) {
      ASSERT_TRUE(ug.degree(v) >= (size_t)M);
    }
  }
}

TEST(gen_hk, random_trial) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = agl::random(1000) + 3;
    V N = M + agl::random(1000);
    double P = agl::random(1000) / 1000.0;
    auto es = generate_hk(N, M, P);

    // Number of edges
    size_t expected_edge_num = M * (M - 1) / 2 + (N - M) * M;
    ASSERT_EQ(es.size(), expected_edge_num);

    G g(es);
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));

    // Check degree
    G ug(make_undirected(es));
    for (V v : ug.vertices()) {
      ASSERT_TRUE(ug.degree(v) >= (size_t)M);
    }
  }
}

TEST(gen_hk, small_case) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = 3;
    V N = M + agl::random(2) + 1;
    double P = agl::random(1000) / 1000.0;
    auto es = generate_hk(N, M, P);

    // Number of edges
    size_t expected_edge_num = M * (M - 1) / 2 + (N - M) * M;
    ASSERT_EQ(es.size(), expected_edge_num);

    G g(es);
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));

    // Check degree
    G ug(make_undirected(es));
    for (V v : ug.vertices()) {
      ASSERT_TRUE(ug.degree(v) >= (size_t)M);
    }
  }
}

TEST(gen_ws, random_trial) {
  uniform_real_distribution<double> rng(1.0);
  for (int trial = 0; trial < 10; ++trial) {
    V N = agl::random(1000) + 10;
    V avg_deg = agl::random(N / 2) + 2;
    if (avg_deg % 2 == 1) avg_deg--;
    double P = rng(agl::random);

    auto es = generate_ws(N, avg_deg, P);
    G g(es);
    pretty_print(g);
  }
}

TEST(gen_config, random_trial) {
  for (int trial = 0; trial < 10; ++trial) {
    V N = agl::random(1000) + 10;
    vector<size_t> deg_seq(N);
    int d = (agl::random(N / 4) + 1) * 2;
    for (int i = 0; i < d * N; ++i) {
      V v;
      do {
        v = agl::random(N);
      } while ((int)deg_seq[v] >= N);
      deg_seq[v]++;
    }

    auto es = generate_config(N, deg_seq);
    G g(es);
    pretty_print(g);

    for (V i = 0; i < N; ++i) {
      ASSERT_TRUE(g.degree(i, kFwd) + g.degree(i, kBwd) <= deg_seq[i]);
    }
  }
}

TEST(gen_kronecker, random_trial) {
  uniform_real_distribution<double> rng(0.4, 0.8);
  for (int trial = 0; trial < 10; ++trial) {
    int scale = agl::random(5) + 6;
    vector<vector<double>> mat(2, vector<double>(2));
    double sum = 0;
    for (int i : make_irange(2)) {
      for (int j : make_irange(2)) {
        mat[i][j] = rng(agl::random);
        sum += mat[i][j];
      }
    }
    G g(generate_kronecker(scale, mat));
    pretty_print(g);

    V num_v = 1 << scale;

    for (int i : make_irange(2)) {
      for (int j : make_irange(2)) {
        mat[i][j] /= sum;
      }
    }

    size_t avg_deg = 16;
    G h(generate_kronecker(scale, avg_deg, mat));
    pretty_print(h);
    ASSERT_TRUE(h.num_edges() <= avg_deg * num_v);
  }
}

TEST(gen_uv_flower, random_trial) {
  for (int trial = 0; trial < 10; ++trial) {
    V u = agl::random(3) + 1;
    V v = u + agl::random(3) + 1;
    V w = u + v;
    V req = w;
    int n = agl::random(3);
    for (int i = 0; i < n; ++i) {
      req *= w;
    }
    cerr << req << " " << u << " " << v << endl;
    auto es = generate_uv_flower(req, u, v);

    // Number of edges
    size_t expected_edge_num = w;
    size_t expected_node_num = w;
    size_t max_deg = 2;
    while (expected_edge_num < es.size()) {
      expected_edge_num *= w;
      expected_node_num = w * expected_node_num - w;
      max_deg *= 2;
    }
    ASSERT_EQ(es.size(), expected_edge_num);

    G g(make_undirected(es));
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));
    ASSERT_EQ(g.num_vertices(), expected_node_num);

    vector<V> deg_check(expected_node_num, 0);
    V current = w;
    for (int i = 0; i < current; ++i) deg_check[i] = 2;
    while (current < expected_node_num) {
      current = current * w - w;
      for (int i = 0; i < current; ++i) {
        if (deg_check[i] == 0)
          deg_check[i] = 2;
        else
          deg_check[i] *= 2;
      }
    }

    // Check degree
    for (int i = 0; i < g.num_vertices(); ++i) {
      ASSERT_EQ(g.degree(i), deg_check[i]);
    }
  }
}

TEST(gen_shm, small_case) {
  V initial_num = 5;
  V required_num = initial_num;
  int t = 2;
  int generation = 3;
  size_t max_deg = initial_num - 1;
  for (int i = 1; i < generation; ++i) {
    required_num = (2 * t + 1) * required_num - 2 * t;
    max_deg *= t;
  }
  auto es = generate_shm(required_num, initial_num, t);

  // Number of edges
  ASSERT_EQ(es.size(), required_num - 1);

  G g(make_undirected(es));
  pretty_print(g);
  ASSERT_TRUE(is_connected(g));
  ASSERT_EQ(g.num_vertices(), required_num);
  ASSERT_EQ(g.degree(0), max_deg);

  // Check degree
  for (V v : g.vertices()) {
    V deg = g.degree(v);
    while (deg % t == 0) deg /= t;
    ASSERT_TRUE(deg == initial_num - 1 || deg == 1 || deg == 2);
  }
}

TEST(gen_shm, random_trial) {
  for (int trial = 0; trial < 10; ++trial) {
    V initial_num = agl::random(20) + 3;
    V required_num = initial_num;
    int t = agl::random(20) + 2;
    int generation = agl::random(5) + 1;
    size_t max_deg = initial_num - 1;
    for (int i = 1; i < generation; ++i) {
      required_num = (2 * t + 1) * required_num - 2 * t;
      max_deg *= t;
    }
    auto es = generate_shm(required_num, initial_num, t);

    // Number of edges
    ASSERT_EQ(es.size(), required_num - 1);

    G g(make_undirected(es));
    pretty_print(g);
    ASSERT_TRUE(is_connected(g));
    ASSERT_EQ(g.num_vertices(), required_num);
    ASSERT_EQ(g.degree(0), max_deg);
    // Check degree
    for (V v : g.vertices()) {
      if (v == 0) continue;
      V deg = g.degree(v);
      while (deg % t == 0) deg /= t;
      ASSERT_TRUE(deg == initial_num - 1 || deg == 1 || deg == 2);
    }
  }
}
