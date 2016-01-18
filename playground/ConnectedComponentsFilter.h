#pragma once

template<class handler_t>
class ConnectedComponentsFilter {
public:
  ConnectedComponentsFilter(const G& g) 
  : n(g.num_vertices()), uf(n), local_indices(n), handlers_indices(n), num_connected_components_(0) {

    FOR(v, n) for (auto e : g.edges(v)) {
      V u = to(e);
      uf.unite(u, v);
    }
    num_connected_components_ = 0;
    FOR(v, n) {
      if (uf.root(v) != v) local_indices[v] = ++local_indices[uf.root(v)];
      else handlers_indices[v] = num_connected_components_++;
    }
    fprintf(stderr, "num_connected_components = %d\n", num_connected_components_);

    vector<bool> used(n);
    FOR(v, n) {
      if (uf.root(v) != v) continue;
      used[v] = true;

      const int num_vs = local_indices[uf.root(v)] + 1;
      local_indices[uf.root(v)] = 0;
      vector<pair<V, V>> edges;
      queue<int> q;
      q.push(v);
      while (!q.empty()) {
        V u = q.front(); q.pop();
        FOR(dir, 2) for (auto& e : g.edges(u, D(dir))) {
          V w = to(e);
          if (!used[w]) {
            used[w] = true;
            q.push(w);
          }
          if (dir == 0) {
            edges.emplace_back(local_indices[u], local_indices[w]);
          }
        }
      }

      fprintf(stderr, "root = %d num_vs = %d, edge_size = %d\n", v, num_vs, sz(edges));
      handlers_.emplace_back(edges, num_vs);
    }
  }

  int query(V u, V v) {
    if (!uf.is_same(u, v)) return 0;
    int lu = local_indices[u], lv = local_indices[v];
    auto& handler = handlers_[handlers_indices[uf.root(u)]];
    if (lu == lv) {
      printf("u = %d, lu = %d, v = %d, lv = %d\n", u, lu, v, lv);
      CHECK(false);
    }
    return handler.query(lu, lv);
  }

  int num_connected_components() const {
    return num_connected_components_;
  }

  const vector<handler_t>& handlers() const {
    return handlers_;
  }

private:
  const int n;
  union_find uf;
  vector<int> local_indices, handlers_indices;
  vector<handler_t> handlers_;
  int num_connected_components_;
};
