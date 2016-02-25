#pragma once
#include "ConnectedComponentsFilter.h"

DEFINE_int32(adjacents_max_deg, 10, "");

class OptimizedGusfieldWith2ECC2 {
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

  //void gusfield_choice_stpair(vector<int>& mincut_order, const vector<int>& degree) {
  //  if(FLAGS_gusfield_choice_stpair_strategy == "sequential"){
  //    return ;
  //  } else if(FLAGS_gusfield_choice_stpair_strategy == "sort_by_degree_ascending") {
  //      sort(mincut_order.begin(),mincut_order.end(),[&degree](const int l,const int r) {
  //        return degree[l] < degree[r];
  //      });
  //  } else if(FLAGS_gusfield_choice_stpair_strategy == "sort_by_degree_desending") {
  //      sort(mincut_order.begin(),mincut_order.end(),[&degree](const int l,const int r) {
  //        return degree[l] > degree[r];
  //      });
  //  } else if(FLAGS_gusfield_choice_stpair_strategy == "random") {
  //    random_shuffle(mincut_order.begin(),mincut_order.end());
  //  } else {
  //    fprintf(stderr, "unrecognized option '-gusfield_choice_stpair_strategy=%s'\n", FLAGS_gusfield_choice_stpair_strategy.c_str());
  //    exit(-1);
  //  }
  //}

  //同じ部分グラフに所属していれば true,そのグループの親nodeを返す
  pair<bool ,V> belong_to_same_component(V s,V t) const {
    if (parent_weight_[t].first == s) swap(s, t); // 下の処理とまとめる
    if (parent_weight_[s].first == t) {
      if (parent_weight_[s].second == -1)
        return make_pair(true, t); // t group
      else return make_pair(false, -1); // 既にコスト確定 = 別のグループである
    }
    if (parent_weight_[t].first == parent_weight_[s].first) {
      return make_pair(true, parent_weight_[t].first);
    } else {
      return make_pair(false, -1);
    }
  }

public:

  OptimizedGusfieldWith2ECC2(vector<pair<V, V>> edges, int num_vs) :
    num_vertices_(num_vs),
    parent_weight_(num_vs, make_pair(-1, -1)),
    depth_(num_vs, -1) {

    vector<int> degree(num_vertices_);
    for (auto& e : edges) degree[e.first]++, degree[e.second]++;

    root_node_ = -1;
    FOR(v, num_vertices_) if (degree[v] != 2) {
      root_node_ = v; break;
    }
    if (root_node_ == -1) root_node_ = 0;

  //phase1
  // size2のカットを見つける
  FOR(v, num_vertices_) {
      if (v == root_node_) continue;
      parent_weight_[v].first = root_node_;
      if (degree[v] == 2) {
        parent_weight_[v].second = 2;
      }
    }
    parent_weight_[root_node_].first = -1;

    dinic_twosided dc_base(edges, num_vs);

  vector<int> used(num_vertices_);
  int used_revision = 0;
  queue<int> q;

  //phase2
  //隣接頂点同士を見る
  for (const auto& e : edges) {
    V s, t; tie(s, t) = e;
    bool same_component; V par; tie(same_component, par) = belong_to_same_component(s, t);
    if (!same_component) continue;

    int min_deg = min(degree[s], degree[t]);
    if(min_deg > FLAGS_adjacents_max_deg) continue;

    dc_base.reset_graph();
    const int cost = dc_base.max_flow(s, t);

    //par -> tへのパスが存在しない => parent nodeがs側にある
    bool parent_belongs_to_s_side = dc_base.path_dont_exists_to_t(par);
    const int F = ++used_revision;

    if (!parent_belongs_to_s_side) {
      // on gomory-fu tree, s -> t -> par
      parent_weight_[s] = make_pair(t, cost);
      //s側に属する頂点の親を,sに変更する
      q.push(s);
      used[s] = F;
      while (!q.empty()) {
        V v = q.front(); q.pop();
        for (auto& e : dc_base.e[v]) {
          const int cap = e.cap(dc_base.graph_revision);
          if (cap == 0 || used[e.to] == F) continue;
          used[e.to] = F;
          q.push(e.to);
          if (parent_weight_[e.to].first == par) parent_weight_[e.to].first = s;
        }
      }
    } else {
      // on gomory-fu tree, t -> s -> par
      parent_weight_[t] = make_pair(s, cost);
      //t側に属する頂点の親を,tに変更する
      q.push(t);
      used[t] = F;
      while (!q.empty()) {
        V v = q.front(); q.pop();
        for (auto& e : dc_base.e[v]) {
          const int cap = dc_base.e[e.to][e.reverse].cap(dc_base.graph_revision);
          if (cap == 0 || used[e.to] == F) continue;
          used[e.to] = F;
          q.push(e.to);
          if (parent_weight_[e.to].first == par) parent_weight_[e.to].first = t;
        }
      }
    }
  }

  //phase3
  //残った頂点をgusfield謹製のアルゴリズムで叩く
  vector<int> mincut_order;
  FOR(s, num_vs) {
    if (s != root_node_ && parent_weight_[s].second == -1) mincut_order.push_back(s);
  }

  for (V s : mincut_order) {
      V t = parent_weight_[s].first;

      dc_base.reset_graph();
      int cost = dc_base.max_flow(s, t);
      // fprintf(stderr, "(%d,%d) cost = %d\n", s, t, cost);
      parent_weight_[s].second = cost;

      if (degree[s] == cost) continue;

      q.push(s);
    const int F = ++used_revision;
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
