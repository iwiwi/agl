#pragma once
#include "ConnectedComponentsFilter.h"
#include "greedy_treepacking.h"

DEFINE_int32(try_greedy_tree_packing, 2, "");
DEFINE_int32(try_large_degree_pairs, 10, "");
DEFINE_int32(choose_near_verices_d, 1, "");
DEFINE_bool(enable_greedy_tree_packing, true, "");
DEFINE_bool(enable_logging_max_flow_details, false, "");
DEFINE_bool(enable_adjacent_cut, true, "");
DEFINE_bool(enable_special_bfs, true, "");

class disjoint_cut_set {
  struct Node {
    int pv, nt;
    int root;
  };

  void erase(int node_id) {
    int pv = nodes[node_id].pv, nt = nodes[node_id].nt;
    if (pv != -1) {
      nodes[pv].nt = nt;
    }
    if (nt != -1) {
      nodes[nt].pv = pv;
    }
    if (pv == -1) {
      root[nodes[node_id].root] = nt;
    }
  }

  void add(int node_id, int group_id) {
    int nt = root[group_id];
    nodes[node_id].root = group_id;
    nodes[node_id].pv = -1;
    nodes[node_id].nt = nt;
    root[group_id] = node_id;
    if (nt != -1) {
      nodes[nt].pv = node_id;
    }
  }

public:
  disjoint_cut_set(int n) : root(n, -1), nodes(n), group_num(1) {
    root[0] = 0;
    nodes[0].pv = -1;
    nodes[n - 1].nt = -1;
    FOR(i, n - 1) {
      nodes[i].nt = i + 1;
      nodes[i + 1].pv = i;
    }
    FOR(i, n) nodes[i].root = 0;
  }

  void create_new_group(int id) {
    erase(id);
    add(id, group_num++);
  }

  bool is_same_group(int a, int b) const {
    if (a >= sz(nodes) || b >= sz(nodes)) return false;
    return nodes[a].root == nodes[b].root;
  }

  void move_other_group(int src, int dst) {
    erase(src);
    add(src, nodes[dst].root);
  }

  bool has_another_id_in_same_group(int id) const {
    const int grp_id = nodes[id].root;
    const int rt = root[grp_id];
    CHECK(rt != -1);
    if (rt != id) return true;
    const int nxt = nodes[rt].nt;
    return nxt != -1;
  }

  int other_id_in_same_group(int id) const {
    const int grp_id = nodes[id].root;
    const int rt = root[grp_id];
    CHECK(rt != -1);
    if (rt != id) return rt;
    const int nxt = nodes[rt].nt;
    CHECK(nxt != -1);
    return nxt;
  }

  pair<int, int> get_two_elements(int group_id) const {
    const int rt = root[group_id];
    CHECK(rt != -1);
    const int nxt = nodes[rt].nt;
    if (nxt == -1) return make_pair(-1, -1);
    return make_pair(rt, nxt);
  }

  bool has_two_elements(int group_id) const {
    auto uv = get_two_elements(group_id);
    return uv.first != -1;
  }

  int debug_group_id(int id) const {
    return nodes[id].root;
  }

  int debug_size_of_group(int grp_id) const {
    int cnt = 0;
    int cur = root[grp_id];
    while (cur != -1) {
      cnt++;
      cur = nodes[cur].nt;
    }
    return cnt;
  }

  void debug() {
    FOR(i, sz(nodes)) {
      int cur = root[i];
      if (cur == -1) continue;
      int cnt = 0;
      while (cur != -1) {
        if (cur >= sz(root)) {
          printf("i = %d, cur = %d\n", i, cur);
        }
        CHECK(cur < sz(root));
        cnt++;
        cur = nodes[cur].nt;
      }
      if (cnt != 1) {
        cur = root[i];
        printf("i = %d : ", i);
        while (cur != -1) {
          printf("%d, ", cur);
          cur = nodes[cur].nt;
        }
        puts("");
      }
    }
  }

  int debug_group_num() const {
    return group_num;
  }

private:
  vector<int> root;
  vector<Node> nodes;
  int group_num;
};

class gomory_hu_tree_builder {
  void dfs(V v, V par = -1) {
    int dep = (par == -1) ? 0 : depth_[par] + 1;
    depth_[v] = dep;
    for (auto& to : edges_[v]) {
      if (to.first == par) continue;
      parent_cost_[to.first] = make_pair(v, to.second);
      dfs(to.first, v);
    }
  }

public:
  gomory_hu_tree_builder(int n) : n_(n), edges_(n), depth_(n), parent_cost_(n) {
    add_edge_count_ = 0;
  }

  void add_edge(V u, V v, int cost) {
    add_edge_count_++;
    edges_[u].emplace_back(v, cost);
    edges_[v].emplace_back(u, cost);
  }

  void build() {
    CHECK(add_edge_count_ == n_ - 1);
    parent_cost_[0] = make_pair(-1, -1);
    dfs(0);
    edges_.clear(); edges_.shrink_to_fit();
  }

  int query(V u, V v) const {
    CHECK(u != v);
    CHECK(u < n_ && v < n_);
    int ans = numeric_limits<int>::max();
    while (u != v) {
      if (depth_[u] > depth_[v]) {
        ans = min(ans, parent_cost_[u].second);
        u = parent_cost_[u].first;
      } else {
        ans = min(ans, parent_cost_[v].second);
        v = parent_cost_[v].first;
      }
    }
    return ans;
  }

  const vector<pair<V, int>>& parent_weight() const {
    return parent_cost_;
  }

  int debug_add_edge_count() const {
    return add_edge_count_;
  }

private:
  int add_edge_count_;
  int n_;
  vector<vector<pair<V, int>>> edges_;
  vector<int> depth_;
  vector<pair<V, int>> parent_cost_;
};

class OptimizedGusfieldWith2ECC {
  void prune_obvious_mincut(vector<pair<V,V>>& edges, disjoint_cut_set& dcs, const vector<int>& degree) {
    vector<int> current_parent(num_vertices_, -1);
    vector<int> current_weight(num_vertices_, -1);
    greedy_treepacking packing_base(edges, num_vertices_);

    auto set_solved = [&](V v, V parent, int weight) {
      current_parent[v] = parent;
      current_weight[v] = weight;
    };

    //degreeの最も大きな頂点をrootに
    int temp_root = 0;
    FOR(v, num_vertices_) {
      if (degree[temp_root] < degree[v]) {
        temp_root = v;
      }
    }

    //次数2のcutを設定
    FOR(v, num_vertices_) {
      if (v == temp_root) continue;
      if (degree[v] == 2) set_solved(v, temp_root, 2); //二重連結成分分解後なので自明なcut
    }

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
        if (degree[v] == 2) continue; // 自明なcutがある

        //debug infomation
        auto gtp_edge_count_before = gtp_edge_count;
        auto gtp_edge_miss_before = gtp_edge_miss;
        auto gtp_edge_use_before = gtp_edge_use;

        auto packing = packing_base;
        packing.arborescence_packing(v);

        //debug infomation
        if (num_vertices_ > 100) { 
          JLOG_ADD("prune_obvious_mincut.gtp_edge_count", gtp_edge_count - gtp_edge_count_before);
          JLOG_ADD("prune_obvious_mincut.gtp_edge_miss", gtp_edge_miss - gtp_edge_miss_before);
          JLOG_ADD("prune_obvious_mincut.gtp_edge_use", gtp_edge_use - gtp_edge_use_before);
        }


        if (current_parent[v] != -1) {
          current_parent[v] = v; // 閉路が出来上がるのを防ぐために、親を自分自身であると登録しておく
        }
        FOR(to, num_vertices_) {
          if (to == v) continue;
          if (current_parent[to] != -1) continue;
          //tree packingの結果がdegreeと一致するなら、flowは流さなくてよい
          if (packing.inedge_count(to) == degree[to]) {
            set_solved(to, v, packing.inedge_count(to));
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

    // cutの求まっていない頂点達について、gusfieldでcutを求める
    int pruned = 0;
    FOR(v, num_vertices_) {
      if (v == temp_root) continue;
      if (current_parent[v] != -1) {
        // cutがもとまっている
        dcs.create_new_group(v);
        gh_builder_.add_edge(v, current_parent[v], current_weight[v]);
        pruned++;
      }
    }

    if (num_vertices_ > 100) {
      JLOG_OPEN("prune") {
        JLOG_ADD("num_vs", num_vertices_);
        JLOG_ADD("pruned", pruned);
      }
    }
  }

  void erase_deg2_edges(vector<pair<V,V>>& edges, vector<int>& degree) {
    const int n = sz(degree);
    vector<vector<int>> e(n);

    for (auto& uv : edges) {
      int u, v; tie(u, v) = uv;
      e[u].push_back(v);
      e[v].push_back(u);
    }

    FOR(i, n) if (sz(e[i]) == 2) {
      int a = e[i][0], b = e[i][1];
      e[a].erase(find(e[a].begin(), e[a].end(), i));
      e[b].erase(find(e[b].begin(), e[b].end(), i));
      e[a].push_back(b);
      e[b].push_back(a);
      e[i].clear();
    }

    edges.clear();

    FOR(i, n) {
      for (auto to : e[i]) {
        if (i < to) {
          edges.emplace_back(i, to);
        }
      }
    }
  }

  queue<int> q;
  vector<int> used;
  int max_flow_times;
  int max_cut_size;
  map<int, int> cutsize_count;
  map<int, int> debug_cut_details;
  vector<int> gomory_hu_cut_used;
  int sep_count;

  void mincut_init() {
    used.resize(num_vertices_);
    max_flow_times = 0;
    max_cut_size = 0;
    while (!q.empty()) q.pop();
    cutsize_count.clear();
    gomory_hu_cut_used.resize(num_vertices_);
    sep_count = 0;
    debug_cut_details.clear();
  }

  void mincut(V s, V t, dinic_twosided& dc_base, disjoint_cut_set& dcs, bool enable_separate_graph = true) {
    if (sz(dc_base.e[s]) > sz(dc_base.e[t])) swap(s, t);

    /*
    * max-flow phase
    */
    const long long before_max_flow = getcap_counter;
    int cost = dc_base.max_flow(s, t);
    const long long after_max_flow = getcap_counter;

    gh_builder_.add_edge(s, t, cost); //cutした結果をgomory_hu treeの枝を登録
                      // fprintf(stderr, "(%d,%d) : %d\n", s, t, cost);
                      //debug infomation
    max_flow_times++;
    if (FLAGS_enable_logging_max_flow_details) {
      JLOG_ADD_OPEN("gusfield.max_flow_details") {
        JLOG_PUT("cost", cost, false);
        JLOG_PUT("edge_count", after_max_flow - before_max_flow, false);
      }
    }
    if (max_flow_times % 10000 == 0) {
      stringstream ss;
      ss << "max_flow_times = " << max_flow_times << ", (" << s << "," << t << ") cost = " << cost;
      JLOG_ADD("gusfield.progress", ss.str());
      fprintf(stderr, "getcap_counter = %lld\n", getcap_counter);
      fprintf(stderr, "addcap_counter = %lld\n", addcap_counter);
      fprintf(stderr, "preflow_eq_degree = %d\n", preflow_eq_degree);
      fprintf(stderr, "flow_eq_0 = %d\n", flow_eq_0);

      fprintf(stderr, "cut details : ");
      for(auto& kv : debug_cut_details) fprintf(stderr, "(%d,%d), ", kv.first, kv.second);
      fprintf(stderr, "\n");
      debug_cut_details.clear();
    }

    //s側の頂点とt側の頂点に分類する
    const int F = max_flow_times;
    int sside = 0, tside = 0;
    if (dc_base.reason_for_finishing_bfs == dinic_twosided::kQsIsEmpty) {
      //s側に属する頂点の親を新しいgroupに移動する
      q.push(s);
      used[s] = F;
      dcs.create_new_group(s);
      while (!q.empty()) {
        V v = q.front(); q.pop();
        sside++;
        for (auto& e : dc_base.e[v]) {
          const int cap = e.cap(dc_base.graph_revision);
          if (cap == 0 || used[e.to] == F) continue;
          used[e.to] = F;
          q.push(e.to);
          if (dcs.is_same_group(t, e.to)) {
            dcs.move_other_group(e.to, s); //tと同じgroupにいた頂点を、s側のgroupに移動
          }
        }
      }
      tside = num_vertices_ - sside;
    } else {
      //t側に属する頂点の親を,tに変更する
      q.push(t);
      used[t] = F;
      dcs.create_new_group(t);
      while (!q.empty()) {
        V v = q.front(); q.pop();
        tside++;
        for (auto& e : dc_base.e[v]) {
          const int cap = dc_base.e[e.to][e.reverse].cap(dc_base.graph_revision);
          if (cap == 0 || used[e.to] == F) continue;
          used[e.to] = F;
          q.push(e.to);
          if (dcs.is_same_group(s, e.to)) {
            dcs.move_other_group(e.to, t); //sと同じgroupにいた頂点を、t側のgroupに移動
          }
        }
      }
      sside = num_vertices_ - tside;
    }

    const int min_side = min(tside, sside);
    if (enable_separate_graph) {
      /*
      * separate phase
      */
      const int separate_upper_bound = 10;
      const bool sside_separate = (sside == min_side && dc_base.reason_for_finishing_bfs == dinic_twosided::kQsIsEmpty && sside > separate_upper_bound);
      const bool tside_separate = (tside == min_side && dc_base.reason_for_finishing_bfs == dinic_twosided::kQtIsEmpty && tside > separate_upper_bound);
      if (sside_separate || tside_separate) {
        sep_count++;
        //gomory_hu algorithm
        //縮約後の頂点2つを追加する
        const int sside_new_vtx = dc_base.n;
        const int tside_new_vtx = sside_new_vtx + 1;
        FOR(_, 2) {
          dc_base.add_vertex();
          used.emplace_back();
          gomory_hu_cut_used.emplace_back();
        }

        int num_reconnected = 0; //枝を繋ぎ直した回数
        if (sside_separate) {
          q.push(s);
          gomory_hu_cut_used[s] = F;
          while (!q.empty()) {
            V v = q.front(); q.pop();
            for (auto& e : dc_base.e[v]) {
              const int cap = e.cap(dc_base.graph_revision);
              if (gomory_hu_cut_used[e.to] == F) continue;
              if (cap == 0) {
                if (used[e.to] != F) {
                  //辺を上手に張り替える
                  dc_base.reconnect_edge(e, sside_new_vtx, tside_new_vtx);
                  num_reconnected++;
                }
              } else {
                gomory_hu_cut_used[e.to] = F;
                q.push(e.to);
              }
            }
          }
        } else {
          q.push(t);
          gomory_hu_cut_used[t] = F;
          while (!q.empty()) {
            V v = q.front(); q.pop();
            for (auto& e : dc_base.e[v]) {
              const int cap = dc_base.e[e.to][e.reverse].cap(dc_base.graph_revision);
              if (gomory_hu_cut_used[e.to] == F) continue;
              if (cap == 0) {
                if (used[e.to] != F) {
                  //辺を上手に張り替える
                  dc_base.reconnect_edge(e, tside_new_vtx, sside_new_vtx);
                  num_reconnected++;
                }
              } else {
                gomory_hu_cut_used[e.to] = F;
                q.push(e.to);
              }
            }
          }
        }
        CHECK(num_reconnected == cost); // 枝を繋ぎ直した回数 == maxflow
      }
    }

    // debug infomation
    cutsize_count[min_side]++;
    debug_cut_details[min_side]++;
    if (min_side > max_cut_size) {
      max_cut_size = min_side;
      if (max_cut_size != 1) {
        fprintf(stderr, "(%d,%d), cut = %d, sside_num = %d, tside_num = %d\n", s, t, cost, sside, tside);
        fprintf(stderr, "  prop : degree[%d] = %d, degree[%d] = %d max_flow_times = %d\n", s, sz(dc_base.e[s]), t, sz(dc_base.e[t]), max_flow_times);
      }
    }
  }

  //次数の大きい頂点対をcutする
  void cut_large_degree_pairs(dinic_twosided& dc_base, disjoint_cut_set& dcs) {
    vector<int> vtxs;
    FOR(v, num_vertices_) {
      if (dcs.has_another_id_in_same_group(v)) vtxs.push_back(v);
    }
    const int tries = max(min(FLAGS_try_large_degree_pairs, sz(vtxs) - 1), 0);
    partial_sort(vtxs.begin(), vtxs.begin() + tries, vtxs.end(), [&dc_base](V l, V r) {
      return sz(dc_base.e[l]) > sz(dc_base.e[r]);
    });

    int cut_large_degree_count = 0;
    for (int i = 1; i <= tries; i++) {
      for (int pari = 0; pari < i; pari++) {
        V s = vtxs[i];
        V t = vtxs[pari];
        if (!dcs.is_same_group(s, t)) continue;
        mincut(s, t, dc_base, dcs);
        cut_large_degree_count++;
      }
    }

    if (cut_large_degree_count > 0) {
      JLOG_PUT("cut_large_degree_pairs.cut_large_degree_count", cut_large_degree_count);
      JLOG_PUT("cut_large_degree_pairs.separated_count", sep_count);
    }
  }

  //隣接頂点同士を見て、まだ切れていなかったらcutする
  void cut_adjacent_pairs(dinic_twosided& dc_base, disjoint_cut_set& dcs) {
    FOR(s, num_vertices_) {
      for(auto& to_edge : dc_base.e[s]) {
        const V t = to_edge.to;
        if (!dcs.is_same_group(s, t)) continue;
        mincut(s, t, dc_base, dcs);
      }
    }

    if (sz(cutsize_count) > 10) {
      stringstream cutsize_count_ss;
      cutsize_count_ss << "cutsize_count : ";
      for (auto& kv : cutsize_count) cutsize_count_ss << "(" << kv.first << "," << kv.second << "), ";
      JLOG_ADD("cut_adjacent_pairs.cutsize_count", cutsize_count_ss.str());
    }
  }

  void gusfield(dinic_twosided& dc_base, disjoint_cut_set& dcs) {
    FOR(group_id, num_vertices_) {
      while (dcs.has_two_elements(group_id)) {
        V s, t; tie(s, t) = dcs.get_two_elements(group_id);
        mincut(s, t, dc_base, dcs);
      }
    }

    if (sz(cutsize_count) > 10) {
      stringstream cutsize_count_ss;
      cutsize_count_ss << "cutsize_count : ";
      for (auto& kv : cutsize_count) cutsize_count_ss << "(" << kv.first << "," << kv.second << "), ";
      JLOG_ADD("gusfield.cutsize_count", cutsize_count_ss.str());
      JLOG_ADD("gomory_hu.separated_count", sep_count);
    }
  }

  void choose_near_verices(dinic_twosided& dc_base, disjoint_cut_set& dcs) {
    vector<int> used;
    used.reserve(num_vertices_ * 2);
    used.resize(dc_base.n, -1);

    FOR(s, num_vertices_) {
      if(!dcs.has_another_id_in_same_group(s)) continue;
      queue<V> q;
      q.push(s);
      used[s] = s;
      FOR(depth, FLAGS_choose_near_verices_d) {
        const int loop_num = sz(q);
        FOR(_, loop_num) {
          const V v = q.front(); q.pop();
          for(auto& to_edge : dc_base.e[v]) {
              const V t = to_edge.to;
              if(used[t] == s) continue;
              used[t] = s;
              q.push(t);
              if(s != t && dcs.is_same_group(s, t)) {
                mincut(s, t, dc_base, dcs);
                used.resize(dc_base.n, -1); // dinic中に頂点数が変わる場合がある
              }
          }
        }
      }
    }

    if (sz(cutsize_count) > 10) {
      stringstream cutsize_count_ss;
      cutsize_count_ss << "cutsize_count : ";
      for (auto& kv : cutsize_count) cutsize_count_ss << "(" << kv.first << "," << kv.second << "), ";
      JLOG_ADD("cut_adjacent_pairs.cutsize_count", cutsize_count_ss.str());
    }
  }

  //次数の最も高い頂点に対して、出来る限りの頂点からflowを流してmincutを求める
  void special_bfs_cut(dinic_twosided& dc_base, disjoint_cut_set& dcs) {
    int max_degree_vtx = 0;
    FOR(v, num_vertices_) if (sz(dc_base.e[max_degree_vtx]) < sz(dc_base.e[v])) max_degree_vtx = v;

    dc_base.special_bfs_init(max_degree_vtx);
    FOR(v, num_vertices_) {
      if (v == max_degree_vtx) continue;
      if (!dcs.is_same_group(v, max_degree_vtx)) continue;
      //graphの形状が変わると損なので、ここでは enable_separate_graph = false する
      mincut(v, max_degree_vtx, dc_base, dcs, false);
    }
  }

public:

  OptimizedGusfieldWith2ECC(vector<pair<V, V>>&& edges, int num_vs) :
    num_vertices_(num_vs),
    gh_builder_(num_vs) {

    if(num_vs > 1000) fprintf(stderr, "OptimizedGusfieldWith2ECC::constructor start : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);
    vector<int> degree(num_vertices_);
    for (auto& e : edges) degree[e.first]++, degree[e.second]++;

    //次数2の頂点と接続を持つ辺を削除して、探索しやすくする
    erase_deg2_edges(edges, degree);

    disjoint_cut_set dcs(num_vs);

    //枝刈り
    if(num_vertices_ > 10000) {
      //頂点数の多いグラフのみlogging
      JLOG_ADD_BENCHMARK("prune_obvious_mincut_time") {
        prune_obvious_mincut(edges, dcs, degree);
      }
    } else {
      prune_obvious_mincut(edges, dcs, degree);
    }

    //debug infomation
    auto preflow_eq_degree_before = preflow_eq_degree;
    auto flow_eq_0_before = flow_eq_0;

    if(num_vs > 1000) fprintf(stderr, "OptimizedGusfieldWith2ECC::dinic_twosided before init : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);
    //dinicの初期化
    dinic_twosided dc_base(std::move(edges), num_vs);
    if(num_vs > 1000) fprintf(stderr, "OptimizedGusfieldWith2ECC::dinic_twosided after init : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);

    mincut_init();

    if (FLAGS_enable_special_bfs) {
      special_bfs_cut(dc_base, dcs);
    }

    // 次数の高い頂点対をcutする
    // グラフをなるべく2分するcutを見つけられると有用
    if (num_vs > 100) {
      cut_large_degree_pairs(dc_base, dcs);
    }

    //まず隣接頂点対からcutしていく
    if (FLAGS_enable_adjacent_cut) {
      CHECK(FLAGS_choose_near_verices_d >= 1);
      if(FLAGS_choose_near_verices_d == 1) {
        cut_adjacent_pairs(dc_base, dcs); 
      } else {
        choose_near_verices(dc_base, dcs);
      }
    }
    // 残った頂点groupをcutする、gomory_hu treeの完成
    gusfield(dc_base, dcs);

    gh_builder_.build();

    auto preflow_eq_degree_after = preflow_eq_degree;
    auto flow_eq_0_after = flow_eq_0;

    if (num_vertices_ > 100) {
      JLOG_OPEN("sp_dfs") {
        JLOG_ADD("preflow_eq_degree", preflow_eq_degree_after - preflow_eq_degree_before);
        JLOG_ADD("flow_eq_0", flow_eq_0_after - flow_eq_0_before);
      }
    }
  }

  int query(V u, V v) const {
    return gh_builder_.query(u, v);
  }

  const vector<pair<V, int>>& parent_weight() const {
    return gh_builder_.parent_weight();
  }

private:
  const int num_vertices_;
  gomory_hu_tree_builder gh_builder_;
};