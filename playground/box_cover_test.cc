#include "box_cover.h"
#include "gtest/gtest.h"
using namespace agl;
using namespace std;

TEST(box_cover, memb) {
  for (int trial = 0; trial < 10; ++trial) {
    V M = 3;
    V N = M + agl::random(1000);
    auto es = generate_ba(N, M);
    G g(es);
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