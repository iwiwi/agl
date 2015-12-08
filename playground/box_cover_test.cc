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
    vector<V> rank(N);
    vector<V> inv(N);
    for (V i = 0; i < N; ++i) {
      inv[i] = i;
    }
    random_shuffle(inv.begin(), inv.end());
    for (int i = 0; i < N; ++i) {
      rank[inv[i]] = i;
    }
    vector<bool> a;
    vector<map<V, V>> naive_x = naive_build_sketch(g, radius, k, rank, a);
    vector<map<V, V>> x = build_sketch(g, radius, k, rank, inv, a);

    for (V v = 0; v < N; v++) {
      ASSERT_EQ(naive_x[v], x[v]) << v;
    }
  }
}
