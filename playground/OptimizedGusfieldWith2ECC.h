#pragma once
#include "ConnectedComponentsFilter.h"
#include "greedy_treepacking.h"

DEFINE_string(gusfield_choice_stpair_strategy, "sequential", "sequential, sort_by_degree_ascending, sort_by_degree_desending, random");
DEFINE_int32(try_greedy_tree_packing, 10, "");
DEFINE_bool(enable_greedy_tree_packing, true, "");

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

  void prune_obvious_mincut(greedy_treepacking& packing, parent_tree_set& pts, vector<int>& mincut_order, const vector<int>& degree) {
    vector<int> current_parent(num_vertices_, -1);

    if (FLAGS_enable_greedy_tree_packing) {
      //次数の高い頂点から順に、 一定回数 greedy tree packingを行って、flowの下界を求める
      const int iteration = min(FLAGS_try_greedy_tree_packing, num_vertices_);
      vector<int> idx(num_vertices_);
      FOR(i, num_vertices_) idx[i] = i;
      partial_sort(idx.begin(), idx.begin() + iteration, idx.end(), [&degree](int l, int r) {
        return degree[l] > degree[r];
      });

      FOR(trying, iteration) {
        const V v = idx[trying];
        packing.arborescence_packing(v);
        if (current_parent[v] != -1) {
          current_parent[v] = v; // 閉路が出来上がるのを防ぐために、親を自分自身であると登録しておく
        }
        FOR(to, num_vertices_) {
          if (to == v) continue;
          if (current_parent[to] != -1) continue;
          //tree packingの結果がdegreeと一致するなら、flowは流さなくてよい
          if (packing.inedge_count(to) == degree[to]) {
            current_parent[to] = v;
            // pts.set_parent(to, v);
            parent_weight_[to].second = packing.inedge_count(to);

            // fprintf(stderr, "(%d, %d), cost = %d\n",v, to, packing.inedge_count(to));
          }
        }
      }

      //閉路が出来上がるのを防ぐためcurrent_parentに代入していた値を、元に戻す
      FOR(trying, iteration) {
        const V v = idx[trying];
        if (current_parent[v] == v) current_parent[v] = -1;
      }

    }

    // root nodeの差し替え
    if (current_parent[root_node_] != -1) {
      pts.change_parent(root_node_, current_parent[root_node_]);
      root_node_ = current_parent[root_node_];
    }

    // cutの求まっていない頂点達について、gusfieldでcutを求める
    FOR(v, num_vertices_) {
      if (v == root_node_) continue;
      if (current_parent[v] != -1) {
        pts.set_parent(v, current_parent[v]);
      } else if (degree[v] == 2) {
        pts.set_parent(v, root_node_);
        parent_weight_[v].second = 2; //ついでに自明なcutを求める
      } else {
        mincut_order.push_back(v);
      }
    }

    if (num_vertices_ > 50) {
      JLOG_OPEN("prune") {
        JLOG_ADD("num_vs", num_vertices_);
        JLOG_ADD("pruned", num_vertices_ - 1 - sz(mincut_order));
      }
    }

  }

public:

  OptimizedGusfieldWith2ECC(vector<pair<V, V>>& edges, int num_vs) :
    num_vertices_(num_vs),
    parent_weight_(num_vs, make_pair(-1, 0)),
    depth_(num_vs, -1) {

    vector<int> degree(num_vertices_);
    for (auto& e : edges) degree[e.first]++, degree[e.second]++;

    //gomory-hu treeのroot nodeの決定
    root_node_ = -1;
    FOR(v, num_vertices_) if (degree[v] != 2) {
      root_node_ = v; break;
    }
    if (root_node_ == -1) root_node_ = 0;
    parent_tree_set pts(num_vs, root_node_);

    vector<int> mincut_order;
    {
      greedy_treepacking packing(edges, num_vs);
      //枝刈りしつつ、まだmincutが求まってない頂点達を見つける
      prune_obvious_mincut(packing, pts, mincut_order, degree);
    }
    //gusfieldのアルゴリズムでcutしていく頂点の順番を決定
    gusfield_choice_stpair(mincut_order, degree);

    //フローを流す gusfieldのアルゴリズム
    dinic_twosided dc_base(edges, num_vs);
    vector<int> used(num_vertices_);
    queue<int> q;
    for (V s : mincut_order) {
      const V t = pts.get_parent(s);

      dc_base.reset_graph();
      int cost = dc_base.max_flow(s, t);
      // fprintf(stderr, "(%d,%d) cost = %d\n", s, t, cost);
      parent_weight_[s].second = cost;

      if (degree[s] == cost) continue;

      //s側のmincutを求めて親の付け替え
      q.push(s);
      const int F = s + 1;
      used[s] = F;
      while (!q.empty()) {
        const V v = q.front(); q.pop();
        for (auto& e : dc_base.e[v]) {
          if (e.cap(dc_base.graph_revision) == 0 || used[e.to] == F) continue;
          used[e.to] = F;
          q.push(e.to);
          const int tpar = pts.get_parent(e.to);
          if (tpar == t) pts.set_parent(e.to, s); //mincut後のe.toはs側に属する
        }
      }
    }

    // gomory-hu treeの親nodeの設定
    FOR(v, num_vs) parent_weight_[v].first = pts.get_parent(v);
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