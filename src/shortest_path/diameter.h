#pragma once

#include <stack>

#include "graph/graph.h"

namespace agl {
template<int kNumDoubleSweep = 10>
int diameter(const G& g) {
  V num_v = g.num_vertices();

  // Decompose the graph into strongly connected components
  std::vector<int> scc(num_v);
  {
    int num_visit = 0, num_scc = 0;
    std::vector<int> ord(num_v, -1);
    std::vector<int> low(num_v);
    std::vector<bool> in(num_v, false);
    std::stack<V> s;
    std::stack<std::pair<V, int>> dfs;

    for (V i : make_irange(num_v)) {
      if (ord[i] != -1) continue;
      dfs.emplace(i, -1);

      while (!dfs.empty()) {
        int v = dfs.top().first;
        int index = dfs.top().second;
        dfs.pop();
        if (index == -1) {
          ord[v] = low[v] = num_visit++;
          s.emplace(v);
          in[v] = true;
        } else {
          low[v] = std::min(low[v], low[to(g.neighbor(v, index))]);
        }

        for (++index; index < g.degree(v); ++index) {
          V w = to(g.neighbor(v, index));
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
  }

  // Compute the diameter lower bound by the double sweep algorithm
  int diameter = 0;
  int qs, qt;
  std::vector<int> dist(num_v, -1);
  std::vector<V> queue(num_v);
  {
    for (int i = 0; i < kNumDoubleSweep; ++i) {
      // forward BFS
      V start = agl::random(num_v);
      qs = qt = 0;
      dist[start] = 0;
      queue[qt++] = start;
      while (qs < qt) {
        V v = queue[qs++];
        for (const auto& e : g.neighbors(v)) {
          if (dist[to(e)] < 0) {
            dist[to(e)] = dist[v] + 1;
            queue[qt++] = to(e);
          }
        }
      }
      for (int j : make_irange(qt)) dist[queue[j]] = -1;

      // barkward BFS
      start = queue[qt - 1];
      qs = qt = 0;
      dist[start] = 0;
      queue[qt++] = start;
      while (qs < qt) {
        V v = queue[qs++];
        for (const auto& e : g.neighbors(v, kBwd)) {
          if (dist[to(e)] < 0) {
            dist[to(e)] = dist[v] + 1;
            queue[qt++] = to(e);
          }
        }
      }
      diameter = std::max(diameter, dist[queue[qt - 1]]);
      for (int j : make_irange(qt)) dist[queue[j]] = -1;
    }
  }

  // Order vertices
  std::vector<std::pair<int64_t, V>> order(num_v);
  {
    for (V v : make_irange(num_v)) {
      int in = 0, out = 0;
      for (const auto& e : g.neighbors(v, kBwd)) {
        if (scc[to(e)] == scc[v]) ++in;
      }
      for (const auto& e : g.neighbors(v)) {
        if (scc[to(e)] == scc[v]) ++out;
      }
      // SCC : reverse topological order
      // inside an SCC : decreasing order of the product of the indegree and outdegree for vertices in the same SCC
      order[v] = std::make_pair(((int64_t)scc[v] << 32) - in * out, v);
    }
    std::sort(order.begin(), order.end());
  }

  // Examine every vertex
  std::vector<int> ecc(num_v, num_v);
  {
    for (V i : make_irange(num_v)) {
      V u = order[i].second;
      if (ecc[u] <= diameter) continue;

      // Refine the eccentricity upper bound
      int ub = 0;
      std::vector<std::pair<int, int>> neighbors;
      for (const auto& e : g.neighbors(u)) neighbors.emplace_back(scc[to(e)], ecc[to(e)] + 1);
      sort(neighbors.begin(), neighbors.end());

      for (size_t j = 0; j < neighbors.size(); ) {
        int component = neighbors[j].first;
        int lb = num_v;
        for (; j < neighbors.size(); ++j) {
          if (neighbors[j].first != component) break;
          lb = std::min(lb, neighbors[j].second);
        }
        ub = std::max(ub, lb);
        if (ub > diameter) break;
      }

      if (ub <= diameter) {
        ecc[u] = ub;
        continue;
      }

      // Conduct a BFS and update bounds
      qs = qt = 0;
      dist[u] = 0;
      queue[qt++] = u;
      while (qs < qt) {
        int v = queue[qs++];
        for (const auto& e : g.neighbors(v)) {
          if (dist[to(e)] < 0) {
            dist[to(e)] = dist[v] + 1;
            queue[qt++] = to(e);
          }
        }
      }
      ecc[u] = dist[queue[qt - 1]];
      diameter = std::max(diameter, ecc[u]);
      for (int j : make_irange(qt)) dist[queue[j]] = -1;

      qs = qt = 0;
      dist[u] = 0;
      queue[qt++] = u;
      while (qs < qt) {
        V v = queue[qs++];
        ecc[v] = std::min(ecc[v], dist[v] + ecc[u]);
        for (const auto& e : g.neighbors(v, kBwd)) {
          // only inside an SCC
          if (dist[to(e)] < 0 && scc[to(e)] == scc[u]) {
            dist[to(e)] = dist[v] + 1;
            queue[qt++] = to(e);
          }
        }
      }
      for (int j : make_irange(qt)) dist[queue[j]] = -1;
    }
  }

  return diameter;
}
}  // namespace agl
