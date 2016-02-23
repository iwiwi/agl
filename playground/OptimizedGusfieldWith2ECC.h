#pragma once
#include "ConnectedComponentsFilter.h"
#include "greedy_treepacking.h"

DEFINE_string(gusfield_choice_stpair_strategy, "sort_by_degree_desending", "sequential, sort_by_degree_ascending, sort_by_degree_desending, random");
DEFINE_int32(try_greedy_tree_packing, 10, "");
DEFINE_bool(enable_greedy_tree_packing, true, "");
DEFINE_bool(enable_logging_max_flow_details, false, "");
DEFINE_bool(enable_adjacent_cut, false, "");

class parent_tree_set {
  int root_;
  vector<int> data_, datainv_, ref_;
  vector<bool> freezed_;
  vector<pair<int, int>> parent_weight_;

public:
  parent_tree_set(const int n, const int root_) : root_(root_), data_(n), datainv_(n), ref_(n, root_),
    freezed_(n), parent_weight_(n, make_pair(-1, -1)) {
    FOR(i, n) datainv_[i] = data_[i] = i;
  }

  int root() const {
    return root_;
  }

  void root(int new_root) {
    root_ = new_root;
  }

  void set_temporary_parent(int v, int par) {
    CHECK(!freezed_[v]);
    ref_[v] = datainv_[par];
  }

  void change_temporary_parent(int cur_parent, int new_parent) {
    int ci = datainv_[cur_parent];
    int ni = datainv_[new_parent];
    swap(data_[ci], data_[ni]);
    swap(datainv_[cur_parent], datainv_[new_parent]);
  }

  int get_temporary_parent(int v) const {
    CHECK(!freezed_[v]);
    return data_[ref_[v]];
  }

  bool is_freezed(int v) const {
    return freezed_[v];
  }

  int get_freezed_parent(int v) const {
    CHECK(freezed_[v]);
    return parent_weight_[v].first;
  }

  void freeze_parent(int v) {
    CHECK(!freezed_[v]);
    parent_weight_[v].first = get_temporary_parent(v);
    freezed_[v] = true;
  }

  void set_weight(int v, int weight) {
    CHECK(!freezed_[v]);
    parent_weight_[v].second = weight;
  }

  // parent_weightをmoveで返して、parent_tree_setを壊します
  vector<pair<int, int>> move_parent_weight() {
    FOR(i, sz(parent_weight_)) if (i != root_) {
      CHECK(freezed_[i]);
    }
    return std::move(parent_weight_);
  }

  void debug() const {
    FOR(v, sz(data_)) {
      printf("%d, ", get_temporary_parent(v));
      if (v != root_ && v == get_temporary_parent(v)) {
        puts("?");
      }
    }
    puts("");
  }
};

class OptimizedGusfieldWith2ECC {
  void build_depth(parent_tree_set &pts) {
    depth_[pts.root()] = 0;

    stack<int> stk;
    FOR(v, num_vertices_) {
      int tmp = v;
      while (depth_[tmp] == -1) {
        stk.push(tmp);
        tmp = parent_weight_[tmp].first;
      }
      while (!stk.empty()) {
        int u = stk.top(); stk.pop();
        depth_[u] = depth_[tmp] + 1;
        tmp = u;
      }
    }
  }

  void gusfield_choice_stpair(parent_tree_set& pts, vector<int>& mincut_order, const vector<int>& degree) {
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

  vector<int> prune_obvious_mincut(parent_tree_set& pts, const vector<int>& degree) {
    vector<int> current_parent(num_vertices_, -1);
    greedy_treepacking packing_base(edges_, num_vertices_);

    auto set_solved = [&](V v, V parent, int weight) {
      current_parent[v] = parent;
      pts.set_weight(v, weight);
    };

    FOR(v, num_vertices_) {
      if (v == pts.root()) continue;
      if (degree[v] == 2) set_solved(v, pts.root(), 2); //二重連結成分分解後なので自明なcut
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
        auto packing = packing_base;
        packing.arborescence_packing(v);
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

    // root nodeの差し替え
    if (current_parent[pts.root()] != -1) {
      int cur_root = pts.root();
      pts.change_temporary_parent(cur_root, current_parent[cur_root]);
      pts.root(current_parent[cur_root]);
    }

    vector<int> mincut_order;
    // cutの求まっていない頂点達について、gusfieldでcutを求める
    FOR(v, num_vertices_) {
      if (v == pts.root()) continue;
      if (current_parent[v] != -1) {
        // cutがもとまっている
        pts.set_temporary_parent(v, current_parent[v]);
        pts.freeze_parent(v);
      } else {
        // もしcutが代入されていた場合に対応するため、-1を代入
        //   「"閉路が出来上がるのを防ぐために、親を自分自身であると登録しておく" の段階で、
        //   既にcutが見つかっていたが、tree packingのroot nodeとなった場合」に、
        //   parentが存在しないのにweightが-1でないという場合が発生する
        pts.set_weight(v, -1);
        mincut_order.push_back(v); //gusfieldで求める
      }
    }

    if (num_vertices_ > 50) {
      JLOG_OPEN("prune") {
        JLOG_ADD("num_vs", num_vertices_);
        JLOG_ADD("pruned", num_vertices_ - 1 - sz(mincut_order));
      }
    }

    return mincut_order;
  }

  void erase_deg2_edges(vector<int>& degree) {
    const int n = sz(degree);
    vector<vector<int>> e(n);

    for (auto& uv : edges_) {
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

    edges_.clear();

    FOR(i, n) {
      for (auto to : e[i]) {
        if (i < to) {
          edges_.emplace_back(i, to);
        }
      }
    }
  }

  //同じ部分グラフに所属していれば true,そのグループの親nodeを返す
  pair<bool, V> belong_to_same_component(parent_tree_set& pts, V s, V t) const {
    if (pts.is_freezed(s) || pts.is_freezed(t)) return make_pair(false, -1); //どちらもコストが決定している

    if (pts.is_freezed(s)) swap(s, t);
    if (pts.is_freezed(t)) { // tはgusfield cut treeの親がfreezeされている
      if (pts.get_temporary_parent(s) == t) {
        return make_pair(true, t); //t == parent
      } else if (pts.get_temporary_parent(s) == pts.get_freezed_parent(t)) {
        return make_pair(true, pts.get_temporary_parent(s));
      } else {
        return make_pair(false, -1); // 同じグループではない
      }
    } else {
      //どちらも親が確定していない
      if (pts.get_temporary_parent(t) == s) swap(s, t);

      if (pts.get_temporary_parent(s) == t) {
        return make_pair(true, t);
      } else if (pts.get_temporary_parent(s) == pts.get_temporary_parent(t)) {
        return make_pair(true, pts.get_temporary_parent(s));
      } else {
        return make_pair(false, -1);
      }
    }
  }


  //隣接頂点同士を見て、まだ切れていなかったらcutする
  void cut_adjacent_pairs(dinic_twosided& dc_base, vector<int>& mincut_order, vector<int>& degree, parent_tree_set& pts) {
    queue<int> q;
    vector<int> used(num_vertices_);
    int used_revision = 0;
    int max_flow_times = 0;
    int mx = 0;
    map<int, int> cutsize_count;

    for (const auto& e : edges_) {
      V s, t; tie(s, t) = e;
      if (degree[s] > degree[t]) swap(s, t);

      bool same_component; V par; tie(same_component, par) = belong_to_same_component(pts, s, t);
      if (!same_component) continue;

      // int min_deg = min(degree[s], degree[t]);
      // if (min_deg > FLAGS_adjacents_max_deg) continue;

      //max-flow
      const long long before_max_flow = getcap_counter;
      dc_base.reset_graph();
      const int cost = dc_base.max_flow(s, t);
      const long long after_max_flow = getcap_counter;

      // fprintf(stderr, "(%d,%d) = %d\n", s, t, cost);

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
      }

      //par -> tへのパスが存在しない => parent nodeがs側にある
      bool parent_belongs_to_s_side = dc_base.path_dont_exists_to_t(par);
      const int F = ++used_revision;

      int sside = 0, tside = 0;
      if (!parent_belongs_to_s_side) {
        // on gomory-fu tree, s -> t -> par
        pts.set_temporary_parent(s, t);
        pts.set_weight(s, cost);
        pts.freeze_parent(s);
        //s側に属する頂点の親を,sに変更する
        q.push(s);
        used[s] = F;
        while (!q.empty()) {
          V v = q.front(); q.pop();
          sside++;
          for (auto& e : dc_base.e[v]) {
            const int cap = e.cap(dc_base.graph_revision);
            if (cap == 0 || used[e.to] == F) continue;
            used[e.to] = F;
            q.push(e.to);
            if (!pts.is_freezed(e.to)) {
              const int tpar = pts.get_temporary_parent(e.to);
              if (tpar == par) pts.set_temporary_parent(e.to, s); //mincut後のe.toはs側に属する
            }
          }
        }
        tside = num_vertices_ - sside;
      } else {
        // on gomory-fu tree, t -> s -> par
        pts.set_temporary_parent(t, s);
        pts.set_weight(t, cost);
        pts.freeze_parent(t);
        //t側に属する頂点の親を,tに変更する
        q.push(t);
        used[t] = F;
        while (!q.empty()) {
          V v = q.front(); q.pop();
          tside++;
          for (auto& e : dc_base.e[v]) {
            const int cap = dc_base.e[e.to][e.reverse].cap(dc_base.graph_revision);
            if (cap == 0 || used[e.to] == F) continue;
            used[e.to] = F;
            q.push(e.to);
            if (!pts.is_freezed(e.to)) {
              const int tpar = pts.get_temporary_parent(e.to);
              if (tpar == par) pts.set_temporary_parent(e.to, t); //mincut後のe.toはt側に属する
            }
          }
        }
        sside = num_vertices_ - tside;
      }

      //debug infomation
      const int x = min(tside, sside);
      cutsize_count[x]++;
      if (x > mx) {
        mx = x;
        if (mx != 1) {
          fprintf(stderr, "(%d,%d), cut = %d, sside_num = %d, tside_num = %d\n", s, t, cost, sside, tside);
          fprintf(stderr, "  prop : degree[%d] = %d, degree[%d] = %d max_flow_times = %d\n", s, degree[s], t, degree[t], max_flow_times);
        }
      }
    }

    if (sz(cutsize_count) > 10) {
      stringstream cutsize_count_ss;
      cutsize_count_ss << "cutsize_count : ";
      for (auto& kv : cutsize_count) cutsize_count_ss << "(" << kv.first << "," << kv.second << "), ";
      JLOG_PUT("cut_adjacent_pairs.cutsize_count", cutsize_count_ss.str());
    }

    //cutしきれなかった頂点たちを、mincut_orderに登録する
    vector<int> mincut_order2;
    for (auto v : mincut_order) {
      if (!pts.is_freezed(v)) mincut_order2.push_back(v);
    }
    mincut_order = mincut_order2;
  }

  void gusfield(dinic_twosided& dc_base, parent_tree_set& pts, vector<int>& mincut_order, vector<int>& degree) {
    vector<int> used(num_vertices_);
    queue<int> q;
    int max_flow_times = 0;
    int mx = 0;

    map<int, int> cutsize_count;
    for (V s : mincut_order) {
      const V t = pts.get_temporary_parent(s);

      //max-flow
      const long long before_max_flow = getcap_counter;
      dc_base.reset_graph();
      int cost = dc_base.max_flow(s, t);
      pts.set_weight(s, cost);
      pts.freeze_parent(s);
      const long long after_max_flow = getcap_counter;

      // fprintf(stderr, "(%d,%d) = %d\n", s, t, cost);

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
      }

      //s側のmincutを求めて親の付け替え
      q.push(s);
      const int F = s + 1;
      used[s] = F;
      int sside = 0;
      while (!q.empty()) {
        const V v = q.front(); q.pop();
        sside++;
        for (auto& e : dc_base.e[v]) {
          if (e.cap(dc_base.graph_revision) == 0 || used[e.to] == F) continue;
          used[e.to] = F;
          q.push(e.to);
          if (!pts.is_freezed(e.to)) {
            const int tpar = pts.get_temporary_parent(e.to);
            if (tpar == t) pts.set_temporary_parent(e.to, s); //mincut後のe.toはs側に属する
          }
        }
      }

      //debug infomation
      int tside = num_vertices_ - sside;
      const int x = min(tside, sside);
      cutsize_count[x]++;
      if (x > mx) {
        mx = x;
        if (mx != 1) {
          fprintf(stderr, "(%d,%d), cut = %d, sside_num = %d, tside_num = %d\n", s, t, cost, sside, tside);
          fprintf(stderr, "  prop : degree[%d] = %d, degree[%d] = %d max_flow_times = %d\n", s, degree[s], t, degree[t], max_flow_times);
        }
      }
    }

    if (sz(cutsize_count) > 10) {
      stringstream cutsize_count_ss;
      cutsize_count_ss << "cutsize_count : ";
      for (auto& kv : cutsize_count) cutsize_count_ss << "(" << kv.first << "," << kv.second << "), ";
      JLOG_PUT("gusfield.cutsize_count", cutsize_count_ss.str());
    }
  }

public:

  OptimizedGusfieldWith2ECC(vector<pair<V, V>>& edges, int num_vs) :
    edges_(edges),
    num_vertices_(num_vs),
    depth_(num_vs, -1) {

    vector<int> degree(num_vertices_);
    for (auto& e : edges) degree[e.first]++, degree[e.second]++;

    //次数2の頂点と接続を持つ辺を削除して、探索しやすくする
    erase_deg2_edges(degree);

    //gomory-hu treeのroot nodeの決定
    int tmp_root_node = -1;
    FOR(v, num_vertices_) if (degree[v] != 2) {
      tmp_root_node = v; break;
    }
    if (tmp_root_node == -1) tmp_root_node = 0;
    parent_tree_set pts(num_vs, tmp_root_node);

    //枝刈りしつつ、まだmincutが求まってない頂点達を見つける
    vector<int> mincut_order = prune_obvious_mincut(pts, degree);

    //dinicの初期化
    dinic_twosided dc_base(edges, num_vs);

    //まず隣接頂点対からcutしていく
    if (FLAGS_enable_adjacent_cut) {
      cut_adjacent_pairs(dc_base, mincut_order, degree, pts);
    }
    //残った頂点gusfieldのアルゴリズムでcutしていく頂点の順番を決定
    gusfield_choice_stpair(pts, mincut_order, degree);
    // gusfieldのアルゴリズムを実行、gomory_hu treeの完成
    gusfield(dc_base, pts, mincut_order, degree);

    // gomory-hu treeの親nodeの設定
    parent_weight_ = pts.move_parent_weight();

    build_depth(pts);
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
  vector<pair<V, V>>& edges_;
  const int num_vertices_;
  vector<pair<V, int>> parent_weight_;
  vector<int> depth_;
};