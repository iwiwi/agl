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
    vector<W> central_distances(g.num_vertices(), g.num_vertices());

    for (V center : memb) {
      vector<bool> vis(g.num_vertices(), false);
      queue<pair<V, W>> que;
      central_distances[center] = 0;
      que.push(make_pair(center, 0));
      while (!que.empty()) {
        V v = que.front().first;
        W dist = que.front().second;
        que.pop();
        for (V u : g.neighbors(v)) {
          if (central_distances[u] <= dist + 1) continue;
          central_distances[u] = dist + 1;
          que.push(make_pair(u, dist + 1));
        }
      }
    }
    for (W dist : central_distances) ASSERT_TRUE(dist <= radius);
  }
}

TEST(box_cover, burning) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = 3;
    V N = M + agl::random(1000);
    auto es = generate_grid(3, 3);
    G g(make_undirected(es));
    pretty_print(g);

    W radius = 1;
    vector<V> burning = box_cover_burning(g, radius);
    vector<W> central_distances(g.num_vertices(), g.num_vertices());

    for (V center : burning) {
      vector<bool> vis(g.num_vertices(), false);
      queue<pair<V, W>> que;
      central_distances[center] = 0;
      que.push(make_pair(center, 0));
      while (!que.empty()) {
        V v = que.front().first;
        W dist = que.front().second;
        que.pop();
        for (V u : g.neighbors(v)) {
          if (central_distances[u] <= dist + 1) continue;
          central_distances[u] = dist + 1;
          que.push(make_pair(u, dist + 1));
        }
      }
    }
    for (W dist : central_distances)
      ASSERT_TRUE(dist <= radius) << dist << " " << radius;
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

    vector<map<V, V>> naive_x = naive_build_sketch(g, radius, k, rank);
    vector<map<V, V>> x = build_sketch(g, radius, k, rank, inv);

    for (V v = 0; v < N; v++) {
      ASSERT_EQ(naive_x[v], x[v]) << v;
    }
  }
}
