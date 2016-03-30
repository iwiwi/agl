#pragma once

DEFINE_int32(goal_oriented_dfs_aster_ub, 2, "");

long long getcap_counter = 0;
long long addcap_counter = 0;
int preflow_eq_degree = 0;
int flow_eq_0 = 0;


class bi_dinitz {
public:
  // two sided bfsが終了した理由
  enum reason_for_finishing_bfs_t {
    kQsIsEmpty,
    kQtIsEmpty,
  };
private:

  class E {
    static const int init_cap_ = 1;
  private: int revision_;
  public:  int to;
  private: int cap_ : 3;
  public:  unsigned int reverse : 29;
  
  public:
    E(int to, int reverse, int cap) :
      revision_(0),to(to),  cap_(cap), reverse(reverse) {
      CHECK(init_cap_ == cap);
    }

    const int cap(int currenct_revision) {
      if (revision_ != currenct_revision) {
        revision_ = currenct_revision;
        cap_ = init_cap_;
      }
      getcap_counter++;
      return cap_;
    }
    void add_cap(int val, int currenct_revision) {
      if (revision_ != currenct_revision) {
        revision_ = currenct_revision;
        cap_ = init_cap_;
      }
      addcap_counter++;
      cap_ += val;
    }

    void reset() {
      revision_ = 0;
      cap_ = init_cap_;
    }
  };

  bool bi_dfs(int s, int t) {
    queue<int> qs, qt;
    qs.push(s); qt.push(t);
    level[s].first = level[t].second = 0;
    bfs_revision[s] = s_side_bfs_revision;
    bfs_revision[t] = t_side_bfs_revision;

    size_t qs_next_get_cap = e[s].size();
    size_t qt_next_get_cap = e[t].size();
    int slevel = 0, tlevel = 0;
    while (qs.size() != 0 && qt.size() != 0) {
      bool path_found = false;
      if (qs_next_get_cap <= qt_next_get_cap) {
        int size = int(qs.size());
        for(int _ = 0; _ < size; _++) {
          const int v = qs.front(); qs.pop();
          qs_next_get_cap -= e[v].size();
          for (auto& t : e[v]) {
            if (t.cap(graph_revision) == 0 || bfs_revision[t.to] == s_side_bfs_revision) continue;
            if (bfs_revision[t.to] == t_side_bfs_revision) {
              path_found = true;
              continue;
            }
            bfs_revision[t.to] = s_side_bfs_revision;
            level[t.to].first = slevel + 1;
            qs_next_get_cap += e[t.to].size();
            qs.push(t.to);
          }
        }
        slevel++;
      } else {
        int size = int(qt.size());
        for(int _ = 0; _ < size; _++) {
          const int v = qt.front(); qt.pop();
          qt_next_get_cap -= e[v].size();
          for (auto& t : e[v]) {
            if (e[t.to][t.reverse].cap(graph_revision) == 0 || bfs_revision[t.to] == t_side_bfs_revision) continue;
            if (bfs_revision[t.to] == s_side_bfs_revision) {
              path_found = true;
              continue;
            }
            bfs_revision[t.to] = t_side_bfs_revision;
            level[t.to].second = tlevel + 1;
            qt_next_get_cap += e[t.to].size();
            qt.push(t.to);
          }
        }
        tlevel++;
      }
      // fprintf(stderr, "slevel : %d, tlevel : %d\n",slevel, tlevel);
      if (path_found) return true;
    }

    reason_for_finishing_bfs = (qs.empty()) ? kQsIsEmpty : kQtIsEmpty;
    return false;
  }

  int dfs(int v, int t, bool use_slevel, int f) {
    // goal_oriented_dfs_aster_ub >= 3 を設定すると、dfs中に同じ辺を使ってしまい、辺のコストが破綻して f < 0 となることがある
    CHECK(f >= 0); 

    if (v == t) return f;
    if (dfs_revision[v] != bfs_revision[v]) {
      dfs_revision[v] = bfs_revision[v];
      iter[v] = 0;
    }
    for (int &i = iter[v]; i < int(e[v].size()); i++) {
      E& _e = e[v][i];
      const int cap = _e.cap(graph_revision);
      if (cap == 0 || bfs_revision[_e.to] / 2 != s_side_bfs_revision / 2) continue;

      bool rec;
      if (use_slevel) rec = bfs_revision[_e.to] == t_side_bfs_revision || level[v].first < level[_e.to].first;
      else rec = bfs_revision[_e.to] == t_side_bfs_revision && level[v].second > level[_e.to].second;
      if (!rec) continue;

      bool next_slevel = use_slevel && bfs_revision[_e.to] == s_side_bfs_revision;
      int d = dfs(_e.to, t, next_slevel, min(f, cap));
      if (d > 0) {
        _e.add_cap(-d, graph_revision);
        e[_e.to][_e.reverse].add_cap(d, graph_revision);
        return d;
      }
    }
    return 0;
  }

  void add_undirected_edge(int f, int t, int c) {
    e[f].push_back(E(t, int(e[t].size()), c));
    e[t].push_back(E(f, int(e[f].size()) - 1, c));
  }

  void reset_rivision() {
    for(int v = 0; v < n; v++) for (auto& e_ : e[v]) e_.reset();
    memset(bfs_revision.data(), 0, sizeof(bfs_revision[0]) * bfs_revision.size());
    memset(dfs_revision.data(), 0, sizeof(dfs_revision[0]) * dfs_revision.size());
    s_side_bfs_revision = 2;
    t_side_bfs_revision = 3;
    graph_revision = 0;
  }

  int goal_oriented_dfs_inner(int v, int flow,int astar_cost) {
    if (v == goal_oriented_bfs_root)
      return flow;

    if (dfs_revision[v] != s_side_bfs_revision) {
      dfs_revision[v] = s_side_bfs_revision;
      iter[v] = 0;
    }
    for (int &i = iter[v]; i < int(e[v].size()); i++) {
      E& to_edge = e[v][i];
      int to = to_edge.to;
      int add_aster_cost = goal_oriented_bfs_depth[to] - goal_oriented_bfs_depth[v] + 1;
      if(add_aster_cost == 2) return 0; //コストの増える頂点は辿らない
      int n_astar_cost = astar_cost + add_aster_cost;
      if (n_astar_cost > FLAGS_goal_oriented_dfs_aster_ub) {
        //sort済なので, これより後で自身の深さよりも浅い頂点は存在しない
        return 0;
      }
      int cap = to_edge.cap(graph_revision);
      if (cap == 0) continue;
      int d = goal_oriented_dfs_inner(to, min(flow, cap), n_astar_cost);
      if (d > 0) {
        to_edge.add_cap(-d, graph_revision);
        e[to_edge.to][to_edge.reverse].add_cap(d, graph_revision);
        return d;
      }
    }

    return 0;
  }

  //v -> goal_oriented_bfs_root にflowを出来る限り送る
  int goal_oriented_dfs(int v) {
    int flow = 0;
    s_side_bfs_revision += 2;
    t_side_bfs_revision += 2;

    for (auto it = e[v].begin(); it != e[v].end(); ++it) {
      auto& to_edge = *it;
      while (to_edge.cap(graph_revision) > 0) {
        int add = goal_oriented_dfs_inner(to_edge.to, to_edge.cap(graph_revision), 0);
        if (add == 0) break;
        flow += add;
        to_edge.add_cap(-add, graph_revision);
        e[to_edge.to][to_edge.reverse].add_cap(add, graph_revision);
      }
    }
    return flow;
  }

public:
  bi_dinitz() : n(0) {}
  bi_dinitz(const G& g)
    : n(g.num_vertices()), level(n), iter(n), bfs_revision(n), dfs_revision(n), e(n),
    s_side_bfs_revision(2), t_side_bfs_revision(3), graph_revision(0), goal_oriented_bfs_root(-1) {
    for(int v = 0; v < n; v++) for (auto& e : g.edges(v)) {
      add_undirected_edge(v, to(e), 1);
    }
  }
  bi_dinitz(const vector<pair<V, V>>& edges, int num_vs)
    : n(num_vs), level(n), iter(n), bfs_revision(n), dfs_revision(n), e(n),
    s_side_bfs_revision(2), t_side_bfs_revision(3), graph_revision(0), goal_oriented_bfs_root(-1) {
    for (auto& uv : edges) {
      add_undirected_edge(uv.first, uv.second, 1);
    }
  }
  bi_dinitz(vector<pair<V, V>>&& edges, int num_vs)
    : n(num_vs), level(n), iter(n), bfs_revision(n), dfs_revision(n), e(n),
    s_side_bfs_revision(2), t_side_bfs_revision(3), graph_revision(0), goal_oriented_bfs_root(-1) {
      edges.shrink_to_fit();

      //こまめに解放しながら辺を追加していく
      while(edges.size() >= 1) {
        size_t loop = max(edges.size() / 2, size_t(10000));
        loop = min(loop, edges.size());
        while(loop--) {
          auto& uv = edges.back();
          add_undirected_edge(uv.first, uv.second, 1);
          edges.pop_back();
        }
        edges.shrink_to_fit();
      }

      for(auto& e_ : e) e_.shrink_to_fit();
  }

  void add_vertex() {
    level.emplace_back();
    iter.emplace_back();
    bfs_revision.emplace_back();
    dfs_revision.emplace_back();
    e.emplace_back();
    n++;
  }

  void reconnect_edge(E& rm, int sside_vtx, int tside_vtx) {
    const int to = rm.to;
    const int to_rev = rm.reverse;
    E& rm_rev = e[to][to_rev];
    const int from = rm_rev.to;
    const int from_rev = rm_rev.reverse;

    add_undirected_edge(sside_vtx, tside_vtx, 1);
    E& se = e[sside_vtx].back();
    E& te = e[tside_vtx].back();

    rm.to = sside_vtx;
    rm.reverse = te.reverse;
    rm_rev.to = tside_vtx;
    rm_rev.reverse = se.reverse;

    se.to = from;
    se.reverse = from_rev;
    te.to = to;
    te.reverse = to_rev;
  }

  int max_flow_core(int s, int t) {
    assert(s != t);
    reset_graph();

    int flow = 0;
    int preflow = 0;
    if (goal_oriented_bfs_root == t) {
      preflow = goal_oriented_dfs(s);
    }

    if(preflow == int(e[s].size())) {
      reason_for_finishing_bfs = kQsIsEmpty;
    } else {
      s_side_bfs_revision += 2;
      t_side_bfs_revision += 2;

      int bfs_counter = 0;
      for (; ; s_side_bfs_revision += 2, t_side_bfs_revision += 2) {
        // fprintf(stderr, "bfs_start\n");
        bool path_found = bi_dfs(s, t);
        bfs_counter++;
        if (!path_found) break;
        while (true) {
          int f = dfs(s, t, true, numeric_limits<int>::max());
          if (f == 0) break;
          flow += f;
        }
      }
      // fprintf(stderr, "bfs_counter : %d\n", bfs_counter);
    }
    // fprintf(stderr, "(%d,%d) : preflow = %d, flow = %d\n", s, t, preflow, flow);
    if(flow == 0 && preflow > 0) {
      if(int(e[s].size()) == preflow) {
        preflow_eq_degree++;
      }
      flow_eq_0++;
    } else {
      // fprintf(stderr, "getcap_counter : %lld\n", getcap_counter);
      // fprintf(stderr, "addcap_counter : %lld\n", addcap_counter);
    }
    return flow + preflow;
  }

  int max_flow(int s, int t) {
    int ans = max_flow_core(s, t);
    return ans;
  }

  bool path_dont_exists_to_t(const int v) const {
    if (reason_for_finishing_bfs == kQsIsEmpty) {
      //sから到達可能な頂点のbfs_revisionには、必ずs_side_bfs_revisionが代入されている
      return bfs_revision[v] == s_side_bfs_revision;
    } else {
      //tから到達不可能
      return bfs_revision[v] != t_side_bfs_revision;
    }
  }

  bool path_dont_exists_from_s(const int v) const {
    if (reason_for_finishing_bfs == kQsIsEmpty) {
      //sから到達不可能
      return bfs_revision[v] != s_side_bfs_revision;
    } else {
      //tから到達可能な頂点のbfs_revisionには、必ずt_side_bfs_revisionが代入されている
      return bfs_revision[v] == t_side_bfs_revision;
    }
  }

  //フローを流す前に実行する
  void reset_graph() {
    graph_revision++;
    if (s_side_bfs_revision >= numeric_limits<decltype(s_side_bfs_revision)>::max() / 2) {
      reset_rivision();
    }
  }

  //rootを起点にbfsをして、 s -> rootのflowを高速化する
  void goal_oriented_bfs_init(int root) {
    goal_oriented_bfs_root = root;
    goal_oriented_bfs_depth.clear();
    goal_oriented_bfs_depth.resize(n, n); // bfsの深さをn(=INF)で初期化

    // set goal_oriented_bfs_depth
    queue<int> q;
    q.push(root);
    goal_oriented_bfs_depth[root] = 0;
    while (!q.empty()) {
      const int v = q.front(); q.pop();
      const int ndepth = goal_oriented_bfs_depth[v] + 1;
      for (auto& to_edge : e[v]) {
        if (goal_oriented_bfs_depth[to_edge.to] > ndepth) {
          goal_oriented_bfs_depth[to_edge.to] = ndepth;
          q.push(to_edge.to);
        }
      }
    }

    //sort edges by depth order
    for(int v = 0; v < n; v++) {
      // dep[l.to]は3つの値しか取らないので高速化可能
      sort(e[v].begin(), e[v].end(), [&dep = goal_oriented_bfs_depth](const E& l, const E& r) {
        return dep[l.to] < dep[r.to];
      });

      //reset reverse edge's "reverse" value
      for(int i = 0; i < int(e[v].size()); i++) {
        const E& to_edge = e[v][i];
        E& rev = e[to_edge.to][to_edge.reverse];
        rev.reverse = i;
      }
    }

  }


  int n;
  vector<pair<int, int>> level;
  vector<int> iter;
  vector<int> bfs_revision, dfs_revision;
  vector<vector<E>> e;
  int s_side_bfs_revision, t_side_bfs_revision;
  int graph_revision;
  reason_for_finishing_bfs_t reason_for_finishing_bfs;

  int goal_oriented_bfs_root;
  vector<int> goal_oriented_bfs_depth;
};
