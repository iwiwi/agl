#pragma once
#include <base/base.h>
#include <graph/graph.h>
#include <vector>
#include <queue>

namespace agl {
template<class handler_t>
class connected_components_filter {
public:
  connected_components_filter(const G& g)
    : n_(g.num_vertices()), uf_(n_), local_indices_(n_), handlers_indices_(n_), num_connected_components_(0) {

    for (int v = 0; v < n_; v++) for (auto e : g.edges(v)) {
      V u = to(e);
      uf_.unite(u, v);
    }
    num_connected_components_ = 0;
    for (int v = 0; v < n_; v++) {
      if (uf_.root(v) != v) local_indices_[v] = ++local_indices_[uf_.root(v)];
      else handlers_indices_[v] = num_connected_components_++;
    }

    std::vector<bool> used(n_);
    for (int v = 0; v < n_; v++) {
      if (uf_.root(v) != v) continue;
      used[v] = true;

      const int num_vs = local_indices_[uf_.root(v)] + 1;
      local_indices_[uf_.root(v)] = 0;
      std::vector<std::pair<V, V>> edges;
      std::queue<int> q;
      q.push(v);
      while (!q.empty()) {
        V u = q.front(); q.pop();
        for (int dir = 0; dir < 2; dir++) for (auto& e : g.edges(u, D(dir))) {
          V w = to(e);
          if (!used[w]) {
            used[w] = true;
            q.push(w);
          }
          if (dir == 0) {
            edges.emplace_back(local_indices_[u], local_indices_[w]);
          }
        }
      }

      edges.shrink_to_fit();
      handlers_.emplace_back(new handler_t(std::move(edges), num_vs));
    }
  }

  int query(V u, V v) {
    if (!uf_.is_same(u, v)) return 0;
    int lu = local_indices_[u], lv = local_indices_[v];
    handler_t& handler = handlers_[handlers_indices_[uf_.root(u)]];
    CHECK(lu != lv);
    return handler.query(lu, lv);
  }

  int num_connected_components() const { return num_connected_components_; }
  const std::vector<std::unique_ptr<handler_t>>& handlers() const { return handlers_; }
  const std::vector<int>& local_indices() const { return local_indices_; }
  const int handlers_index(int v) { return handlers_indices_[uf_.root(v)]; }

private:
  const int n_;
  union_find uf_;
  std::vector<int> local_indices_, handlers_indices_;
  std::vector<std::unique_ptr<handler_t>> handlers_;
  int num_connected_components_;
};
} //namespace agl