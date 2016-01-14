#pragma once
#include "ConnectedComponentsFilter.h"

class OptimizedGusfieldWith2ECC_core {
  void build_depth() {
    depth_[root_node_] = 0;

    stack<int> stk;
    FOR(v, num_vertices_) {
      while (depth_[v] == -1) {
        stk.push(v);
        v = parent_weight_[v].first;
      }
      while (!stk.empty()) {
        int u = stk.top(); stk.pop();
        depth_[u] = depth_[v] + 1;
        v = u;
      }
    }
  }

public:

  OptimizedGusfieldWith2ECC_core(vector<pair<V, V>>& edges, int num_vs) :
    num_vertices_(num_vs),
    parent_weight_(num_vs, make_pair(-1, 0)),
    depth_(num_vs, -1) {

    vector<int> degree(num_vertices_);
    for (auto& e : edges) degree[e.first]++, degree[e.second]++;

    root_node_ = -1;
    FOR(v, num_vertices_) if (degree[v] != 2) {
      root_node_ = v; break;
    }
    if (root_node_ == -1) root_node_ = 0;

    vector<int> mincut_order;
    FOR(v, num_vertices_) {
      if (v == root_node_) continue;
      parent_weight_[v].first = root_node_;
      if (degree[v] == 2) {
        parent_weight_[v].second = 2;
      } else {
        mincut_order.push_back(v);
      }
    }
    parent_weight_[root_node_].first = -1;

    dinic_twosided dc_base(edges, num_vs);

    vector<int> used(num_vertices_);
    for (V s : mincut_order) {
      V t = parent_weight_[s].first;

      auto dc = dc_base;
      int cost = dc.max_flow(s, t);
      fprintf(stderr, "(%d,%d) cost = %d\n", s, t, cost);
      parent_weight_[s].second = cost;

      if (degree[s] == cost) continue;

      queue<int> q;
      q.push(s);
      const int F = s + 1;
      used[s] = F;
      while (!q.empty()) {
        V v = q.front(); q.pop();
        for (auto& e : dc.e[v]) {
          if (e.cap(dc.graph_revision) == 0 || used[e.to] == F) continue;
          used[e.to] = F;
          q.push(e.to);
          if (parent_weight_[e.to].first == t) parent_weight_[e.to].first = s;
        }
      }
    }

    build_depth();
  }

  int query(V u, V v) {
    CHECK(u != v);
    CHECK(u < num_vertices_ && v < num_vertices_);
    int ans = numeric_limits<int>::max();
    while (u != v) {
      if (depth_[u] > depth_[v]) {
        ans = min(ans, parent_weight_[u].second);
        u = parent_weight_[u].first;
      } else {
        ans = min(ans, parent_weight_[v].second);
        v = parent_weight_[v].first;
      }
    }
    return ans;
  }

private:
  const int num_vertices_;
  vector<pair<V, int>> parent_weight_;
  vector<int> depth_;

  int root_node_;
};

// 2ECC = 2-edge connected components
class OptimizedGusfieldWith2ECC {
public:
  const int n;
  const G& g;

  void lowlink_dfs(int v, int par, int& cur_ord) {
    lowlink[v] = order[v] = cur_ord++;
    FOR(dir, 2) for (auto to : g.edges(v, D(dir))) {
      if (to == par) continue;
      if (order[to] == -1) {
        lowlink_dfs(to, v, cur_ord);
        lowlink[v] = min(lowlink[v], lowlink[to]);
        if (order[v] < lowlink[to]) bridge.emplace_back(v, to);
        else biconnected_graphs_edges.emplace_back(v, to);
      } else {
        lowlink[v] = min(lowlink[v], lowlink[to]);
        if (v < to) biconnected_graphs_edges.emplace_back(v, to);
      }
    }
  }

  OptimizedGusfieldWith2ECC(const G& g) : n(g.num_vertices()), g(g), uf(n), lowlink(n, -1), order(n, -1) {
    FOR(v, n) for (auto& e : g.edges(v)) {
      uf.unite(v, to(e));
    }

    const int num_edges = g.num_edges();

    FOR(v, n) if (uf.root(v) == v) {
      int cur_ord = 0;
      lowlink_dfs(v, -1, cur_ord);
    }

    CHECK(sz(bridge) + sz(biconnected_graphs_edges) == num_edges);

    G new_g(biconnected_graphs_edges, n);
    biconnected_graph_handler.reset(new ConnectedComponentsFilter<OptimizedGusfieldWith2ECC_core>(new_g));
  }

public:
  int query(V u, V v) {
    int ans = biconnected_graph_handler->query(u, v);
    if (ans == 0) {
      if (uf.is_same(u, v)) return 1; // 橋で間接的につながっている
      else return 0;
    }
    return ans;
  }

private:
  union_find uf;
  vector<int> lowlink, order;
  vector<pair<V, V>> bridge, biconnected_graphs_edges;

  unique_ptr<ConnectedComponentsFilter<OptimizedGusfieldWith2ECC_core>> biconnected_graph_handler;
};