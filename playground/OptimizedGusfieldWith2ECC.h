#pragma once
#include "ConnectedComponentsFilter.h"

DEFINE_string(gusfield_choice_stpair_strategy, "sequential", "sequential, sort_by_degree_ascending, sort_by_degree_desending, random");

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
    if(FLAGS_gusfield_choice_stpair_strategy == "sequential"){
      return ;
    } else if(FLAGS_gusfield_choice_stpair_strategy == "sort_by_degree_ascending") {
        sort(mincut_order.begin(),mincut_order.end(),[&degree](const int l,const int r) {
          return degree[l] < degree[r];
        });
    } else if(FLAGS_gusfield_choice_stpair_strategy == "sort_by_degree_desending") {
        sort(mincut_order.begin(),mincut_order.end(),[&degree](const int l,const int r) {
          return degree[l] > degree[r];
        });
    } else if(FLAGS_gusfield_choice_stpair_strategy == "random") {
      random_shuffle(mincut_order.begin(),mincut_order.end());
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

    gusfield_choice_stpair(mincut_order, degree);

    dinic_twosided dc_base(edges, num_vs);
    vector<int> used(num_vertices_);
    queue<int> q;
    for (V s : mincut_order) {
      V t = parent_weight_[s].first;

      dc_base.reset_graph();
      int cost = dc_base.max_flow(s, t);
      // fprintf(stderr, "(%d,%d) cost = %d\n", s, t, cost);
      parent_weight_[s].second = cost;

      if (degree[s] == cost) continue;

      q.push(s);
      const int F = s + 1;
      used[s] = F;
      while (!q.empty()) {
        V v = q.front(); q.pop();
        for (auto& e : dc_base.e[v]) {
          if (e.cap(dc_base.graph_revision) == 0 || used[e.to] == F) continue;
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

  const vector<pair<V,int>>& parent_weight() const {
    return parent_weight_;
  }

private:
  const int num_vertices_;
  vector<pair<V, int>> parent_weight_;
  vector<int> depth_;

  int root_node_;
};
