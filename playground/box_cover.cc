#include "box_cover.h"
using namespace std;

vector<V> box_cover_memb(const G &g, W radius) {
  V num_v = g.num_vertices();
  vector<bool> covered(num_v, false);
  vector<bool> centered(num_v, false);
  vector<V> box_centers;

  for (;;) {
    V next = -1, max_mass = 0;
    vector<V> next_nodes;
    for (V pv = 0; pv < num_v; ++pv) {
      if (centered[pv]) continue;

      queue<pair<V, W>> que;
      vector<bool> vis(num_v, false);
      vector<V> nodes;

      que.push(make_pair(pv, 0));
      nodes.push_back(pv);
      vis[pv] = true;

      while (!que.empty()) {
        V v = que.front().first;
        W dist = que.front().second;
        que.pop();
        if (dist >= radius) continue;
        for (V u : g.neighbors(v)) {
          if (vis[u]) continue;
          que.push(make_pair(u, dist + 1));
          nodes.push_back(u);
          vis[u] = true;
        }
      }

      if ((V)nodes.size() > max_mass) {
        next = pv;
        max_mass = nodes.size();
        next_nodes.swap(nodes);
      }
    }
    if (next == -1) break;

    centered[next] = true;
    box_centers.push_back(next);
    for (V node : next_nodes) covered[node] = true;
    if (find(covered.begin(), covered.end(), false) == covered.end()) break;
  }

  return box_centers;
}

vector<V> box_cover_burning(const G &g, W radius) { return {}; }
