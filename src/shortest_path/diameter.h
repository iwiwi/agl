#pragma once

#include <stack>

#include "dijkstra.h"

namespace agl {
template<typename GraphT = G, int kNumDoubleSweep = 10>
typename GraphT::W diameter(const GraphT& g) {
  using W = typename GraphT::W;
  constexpr W INF = std::numeric_limits<W>::max();
  auto cmp = [](W a, W b){return (a != INF ? a : -1) < (b != INF ? b : -1);};
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
        V v = dfs.top().first;
        int index = dfs.top().second;
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
  }

  // Order vertices
  std::vector<std::pair<int64_t, V>> order(num_v);
  {
    for (V v : make_irange(num_v)) {
      int in = 0, out = 0;
      for (const auto& e : g.edges(v, kBwd)) {
        if (scc[to(e)] == scc[v]) ++in;
      }
      for (const auto& e : g.edges(v)) {
        if (scc[to(e)] == scc[v]) ++out;
      }
      // SCC : reverse topological order
      // inside an SCC : decreasing order of the product of the indegree and outdegree for vertices in the same SCC
      order[v] = std::make_pair(((int64_t)scc[v] << 32) - in * out, v);
    }
    std::sort(order.begin(), order.end());
  }

  // Compute the diameter lower bound by the double sweep algorithm
  W diameter = 0;
  {
    for (int i = 0; i < kNumDoubleSweep; ++i) {
      V start = agl::random(num_v);
      auto dist = single_source_distance(g, start);
      start = max_element(dist.begin(), dist.end(), cmp) - dist.begin();
      dist = single_source_distance(g, start);
      diameter = std::max(diameter, *max_element(dist.begin(), dist.end(), cmp));
    }
  }

  // Examine every vertex
  std::vector<W> ecc(num_v, INF);
  {
    for (int i : make_irange(num_v)) {
      V u = order[i].second;
      if (ecc[u] <= diameter) continue;

      // Refine the eccentricity upper bound
      W ub = 0;
      std::vector<std::pair<int, W>> neighbors;
      for (const auto& e : g.edges(u)) {
        if (ecc[to(e)] != INF) {
          neighbors.emplace_back(scc[to(e)], ecc[to(e)] + weight(e));
        }
      }
      sort(neighbors.begin(), neighbors.end());

      for (size_t j = 0; j < neighbors.size(); ) {
        int component = neighbors[j].first;
        W lb = INF;
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

      // Update bounds
      auto dist = single_source_distance(g, u);
      ecc[u] = *max_element(dist.begin(), dist.end(), cmp);
      diameter = std::max(diameter, ecc[u]);

      auto h = make_dijkstra_heap(g);
      h.decrease(u, 0);
      while (!h.empty()) {
        V v = h.top_vertex();
        W w = h.top_weight();
        h.pop();
        ecc[v] = std::min(ecc[v], w + ecc[u]);
        for (const auto &e : g.edges(v, kBwd)) {
          if (scc[to(e)] == scc[u]) {
            h.decrease(to(e), w + weight(e));
          }
        }
      }
    }
  }

  return diameter;
}
}  // namespace agl
