#pragma once
#include "ConnectedComponentsFilter.h"

DEFINE_string(gusfield_choice_stpair_strategy, "sequential", "sequential, sort_by_degree_ascending, sort_by_degree_desending, random");

class parent_tree_set {
  const int root;
  vector<int> data, datainv, ref;

public:
  parent_tree_set(const int n, const int root) : root(root), data(n), datainv(n), ref(n, root) {
    FOR(i, n) datainv[i] = data[i] = i;
  }

  void set_parent(int v, int par) {
    ref[v] = datainv[par];
  }

  void change_parent(int cur_parent, int new_parent) {
    int ci = datainv[cur_parent];
    int ni = datainv[new_parent];
    swap(data[ci], data[ni]);
    swap(datainv[cur_parent], datainv[new_parent]);
  }

  int get_parent(int v) const {
    return data[ref[v]];
  }

  void debug() const {
    FOR(v, sz(data)) {
      printf("%d, ", get_parent(v));
      if (v != root && v == get_parent(v)) {
        puts("?");
      }
    }
    puts("");
  }
};

class OptimizedGusfieldWith2ECC {
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

  void gusfield_choice_stpair(vector<int>& mincut_order, const vector<int>& degree) {
    if (FLAGS_gusfield_choice_stpair_strategy == "sequential") {
      return;
    } else if (FLAGS_gusfield_choice_stpair_strategy == "sort_by_degree_ascending") {
      sort(mincut_order.begin(), mincut_order.end(), [&degree](const int l, const int r) {
        return degree[l] < degree[r];
      });
    } else if (FLAGS_gusfield_choice_stpair_strategy == "sort_by_degree_desending") {
      sort(mincut_order.begin(), mincut_order.end(), [&degree](const int l, const int r) {
        return degree[l] > degree[r];
      });
    } else if (FLAGS_gusfield_choice_stpair_strategy == "random") {
      random_shuffle(mincut_order.begin(), mincut_order.end());
    } else {
      fprintf(stderr, "unrecognized option '-gusfield_choice_stpair_strategy=%s'\n", FLAGS_gusfield_choice_stpair_strategy.c_str());
      exit(-1);
    }
  }


public:

  OptimizedGusfieldWith2ECC(vector<pair<V, V>>& edges, int num_vs) :
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
    parent_tree_set pts(num_vs, root_node_);


    vector<int> mincut_order;
    FOR(v, num_vertices_) {
      if (v == root_node_) continue;
      if (degree[v] == 2) {
        parent_weight_[v].second = 2;
      } else {
        mincut_order.push_back(v);
      }
    }

    gusfield_choice_stpair(mincut_order, degree);

    dinic_twosided dc_base(edges, num_vs);
    vector<int> used(num_vertices_);
    queue<int> q;
    for (V s : mincut_order) {
      V t = pts.get_parent(s);

      dc_base.reset_graph();
      int cost = dc_base.max_flow(s, t);
      // fprintf(stderr, "(%d,%d) cost = %d\n", s, t, cost);
      parent_weight_[s].second = cost;

      if (degree[s] == cost) continue;

      const int F = s + 1;
      if (dc_base.reason_for_finishing_bfs == dinic_twosided::kQsIsEmpty) {
        q.push(s);
        used[s] = F;
        while (!q.empty()) {
          V v = q.front(); q.pop();
          for (auto& e : dc_base.e[v]) {
            if (e.cap(dc_base.graph_revision) == 0 || used[e.to] == F) continue;
            used[e.to] = F;
            q.push(e.to);
            int tpar = pts.get_parent(e.to);
            if (tpar == t) pts.set_parent(e.to, s);
          }
        }
      } else {
        pts.change_parent(s, t);
        pts.set_parent(s, t);
        q.push(t);
        used[t] = F;
        while (!q.empty()) {
          V v = q.front(); q.pop();
          for (auto& e : dc_base.e[v]) {
            const int cap = dc_base.e[e.to][e.reverse].cap(dc_base.graph_revision);
            if (cap == 0 || used[e.to] == F) continue;
            used[e.to] = F;
            q.push(e.to);
            int tpar = pts.get_parent(e.to);
            if (tpar == s) pts.set_parent(e.to, t);
          }
        }
      }

      // pts.debug();
    }

    FOR(v, num_vs) {
      parent_weight_[v].first = pts.get_parent(v);
    }
    parent_weight_[root_node_].first = -1;

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

  const vector<pair<V, int>>& parent_weight() const {
    return parent_weight_;
  }

private:
  const int num_vertices_;
  vector<pair<V, int>> parent_weight_;
  vector<int> depth_;

  int root_node_;
};
