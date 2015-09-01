#pragma once
#include <queue>
#include <algorithm>
#include "graph/graph.h"

namespace agl {
template<typename GraphType>
bool is_connected(const GraphType &g) {
  if (g.num_vertices() == 0) return true;

  std::queue<V> que;
  std::vector<bool> vis(g.num_vertices());
  que.push(0);
  vis[0] = true;
  while (!que.empty()) {
    V v = que.front();
    que.pop();
    for (V tv : g.neighbors(v)) {
      if (vis[tv]) continue;
      que.push(tv);
      vis[tv] = true;
    }
  }
  return std::find(vis.begin(), vis.end(), false) == vis.end();
}
}  // namespace agl
