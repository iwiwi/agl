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
    ASSERT_TRUE(g.num_vertices() == num_v);

    for (int i : make_irange(2)) {
      for (int j : make_irange(2)) {
        mat[i][j] /= sum;
      }
    }

    size_t avg_deg = 16;
    G h(generate_kronecker(scale, avg_deg, mat));
    pretty_print(h);
    ASSERT_TRUE(h.num_vertices() == num_v);
    ASSERT_TRUE(h.num_edges() <=  avg_deg * num_v);
  }
}
