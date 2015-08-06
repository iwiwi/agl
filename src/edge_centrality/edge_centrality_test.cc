#include "edge_centrality.h"
#include "gtest/gtest.h"
using namespace std;
using namespace agl;

TEST(edge_betweenness_centrality, sample_single) {
  FLAGS_edge_betweenness_centrality_num_samples = 1;

  for (int t = 0; t < 10000; ++t) {
    G g(gen_erdos_renyi(8, 3));
    // pretty_print(g);

    const edge_centrality_map out1 = edge_betweenness_centrality_sample_slow(g);
    const edge_centrality_map out2 = edge_betweenness_centrality_sample(g);

    for (V u : g.vertices()) {
      for (E e : g.edges(u)) {
        V v = to(e);
        double o1 = out1.at(make_pair(u, v));
        double o2 = out2.at(make_pair(u, v));
        // printf("%d->%d: %f %f\n", u, v, o1, o2);
        ASSERT_TRUE(is_eq(o1, o2));
      }
    }
  }
}

TEST(edge_betweenness_centrality, sample_multiple) {
  FLAGS_edge_betweenness_centrality_num_samples = 2;

  for (int t = 0; t < 1000; ++t) {
    G g(gen_erdos_renyi(8, 2));
    // pretty_print(g);

    const edge_centrality_map out1 = edge_betweenness_centrality_sample_slow(g);
    const edge_centrality_map out2 = edge_betweenness_centrality_sample(g);

    for (V u : g.vertices()) {
      for (E e : g.edges(u)) {
        V v = to(e);
        double o1 = out1.at(make_pair(u, v));
        double o2 = out2.at(make_pair(u, v));
        // printf("%d->%d: %f %f\n", u, v, o1, o2);
        ASSERT_TRUE(is_eq(o1, o2));
      }
    }
  }
}

TEST(edge_betweenness_centrality, random) {
  constexpr int kNumSamples = 10000;
  constexpr double kErrorTolerance = 5 * 1 / sqrt(kNumSamples);
  FLAGS_edge_betweenness_centrality_num_samples = kNumSamples;

  for (int t = 0; t < 100; ++t) {
    G g(gen_erdos_renyi(30, 2));
    const edge_centrality_map ans = edge_betweenness_centrality_naive(g);
    const edge_centrality_map out1 = edge_betweenness_centrality_sample_slow(g);
    const edge_centrality_map out2 = edge_betweenness_centrality_sample(g);

    for (V u : g.vertices()) {
      for (E e : g.edges(u)) {
        V v = to(e);
        double a = ans.at(make_pair(u, v));
        double o1 = out1.at(make_pair(u, v));
        double o2 = out2.at(make_pair(u, v));

        ASSERT_TRUE(is_eq(o1, o2));
        ASSERT_LT(fabs(a - o1), kErrorTolerance);
        ASSERT_LT(fabs(a - o2), kErrorTolerance);
      }
    }
  }
}
