#pragma once

template<class max_flowT>
class OptimizedGusfield {
  void build_depth() {
    depth_[0] = 0;

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
  OptimizedGusfield(vector<pair<V, V>>& edges, int num_vs) :
    num_vertices_(num_vs),
    parent_weight_(num_vs, make_pair(-1, 0)),
    depth_(num_vs, -1) {
    for (V v = 1; v < num_vertices_; v++) {
      parent_weight_[v].first = 0;
    }
    max_flowT dc_base(edges, num_vs);

    vector<int> degree(num_vertices_);
    for (auto& e : edges) degree[e.first]++, degree[e.second]++;

    vector<int> used(num_vertices_);
    FOR(s, num_vertices_) {
      if (parent_weight_[s].first == -1) continue;
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
};