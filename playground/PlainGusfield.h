#pragma once
#include "dinitz.h"

class PlainGusfield{
  int query_dfs(V v, V t, int cost, V par = -1) const {
    if (v == t) return cost;
    for (const auto& to_cost : binary_tree_edges_[v]) {
      V to; int edge_cost; tie(to, edge_cost) = to_cost;
      if (to == par) continue;
      int ncost = query_dfs(to, t, min(cost, edge_cost), v);
      if (ncost != -1) return ncost;
    }

    return -1;
  }

public:

  void add_edge(V s,V t,int cost){
      binary_tree_edges_[s].emplace_back(t, cost);
      binary_tree_edges_[t].emplace_back(s, cost);
  }

  PlainGusfield(G& g) : num_vertices_(g.num_vertices()), binary_tree_edges_(g.num_vertices()) {
    union_find uf(num_vertices_);

    FOR(v, num_vertices_) for (auto& e : g.edges(v)) {
      uf.unite(v, to(e));
    }
    vector<int> p(num_vertices_);
    FOR(v, num_vertices_) {
      if (v == uf.root(v)) p[v] = -1;
      else p[v] = uf.root(v);
    }
    vector<int> root_vtxs;
    FOR(v, num_vertices_) if(p[v] == -1) root_vtxs.push_back(v);
    FOR(i, sz(root_vtxs) - 1) {
      add_edge(root_vtxs[i], root_vtxs[i+1], 0);
    }

    FOR(s, num_vertices_) {
      if (p[s] == -1) continue;
      V t = p[s];
      dinitz dc(num_vertices_);
      FOR(v, num_vertices_) {
        if (!uf.is_same(s, v)) continue;
        for (auto& e : g.edges(v)) {
          dc.add_undirected_edge(v, to(e), 1);
        }
      }

      int cost = dc.max_flow(s, t);
      // fprintf(stderr, "(%d,%d) cost = %d\n", s, t, cost);
      add_edge(s, t, cost);


      vector<char> used(num_vertices_);
      queue<int> q;
      q.push(s);
      used[s] = true;
      while (!q.empty()) {
        V v = q.front(); q.pop();
        for (auto& e : dc.e[v]) {
          if (e.cap == 0 || used[e.to]) continue;
          used[e.to] = true;
          q.push(e.to);
          if (p[e.to] == t) p[e.to] = s;
        }
      }
    }
  }
  
  int query(V u, V v) const {
    CHECK(u != v);
    CHECK(u < num_vertices_ && v < num_vertices_);
    int ans = query_dfs(u, v, numeric_limits<int>::max());
    if (ans == -1) {
      return 0; // 到達できなかった
    }
    return ans;
  }

  void print_gomory_hu_tree_dfs(V v,V par,ostream& os) {
    for(auto& to : binary_tree_edges_[v]) {
      if(to.first == par) continue;
      os << v << " " << to.first << " " << to.second << endl;
      print_gomory_hu_tree_dfs(to.first, v, os);
    }
  }

  void print_gomory_hu_tree(ostream& os) {
    print_gomory_hu_tree_dfs(0, -1, os);
  }

private:
  const int num_vertices_;
  vector<vector<pair<V, int>>> binary_tree_edges_;
};