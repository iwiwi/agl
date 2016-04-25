#pragma once

#include <stack>

#include "graph/graph.h"

namespace agl {
template<typename GraphT = G>
std::vector<V> strongly_connected_components(const GraphT& g) {
  const V num_v = g.num_vertices();
  std::vector<V> scc(num_v);
  V num_visit = 0, num_scc = 0;
  std::vector<V> ord(num_v, -1);
  std::vector<V> low(num_v);
  std::vector<bool> in(num_v, false);
  std::stack<V, std::vector<V>> s;
  std::stack<std::pair<V, V>, std::vector<std::pair<V, V>>> dfs;

  for (V i : make_irange(num_v)) {
    if (ord[i] != -1) continue;
    dfs.emplace(i, -1);

    while (!dfs.empty()) {
      V v = dfs.top().first;
      V index = dfs.top().second;
      dfs.pop();
      if (index == -1) {
        ord[v] = low[v] = num_visit++;
        s.emplace(v);
        in[v] = true;
      } else {
        low[v] = std::min(low[v], low[to(g.edge(v, index))]);
      }

      for (++index; index < g.degree(v); ++index) {
        V w = to(g.edge(v, index));
        if (ord[w] == -1) {
          dfs.emplace(v, index);
          dfs.emplace(w, -1);
          break;
        } else if (in[w]) {
          low[v] = std::min(low[v], ord[w]);
        }
      }

      if (index == g.degree(v) && low[v] == ord[v]) {
        while (true) {
          V w = s.top();
          s.pop();
          in[w] = false;
          scc[w] = num_scc;
          if (v == w) break;
        }
        ++num_scc;
      }
    }
  }
  return scc;
}
}  // namespace agl
