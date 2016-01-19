#pragma once

class dinic_twosided {
  class E {
    const int init_cap_;
    int cap_, revision_;
  public:
    E(int to, int reverse, int cap) :
      init_cap_(cap), cap_(cap), revision_(0), to(to), reverse(reverse) {
    }

    const int to, reverse;
    const int cap(int currenct_revision) {
      if (revision_ != currenct_revision) {
        revision_ = currenct_revision;
        cap_ = init_cap_;
      }
      return cap_;
    }
    void add_cap(int val, int currenct_revision) {
      revision_ = currenct_revision;
      cap_ += val;
    }

    void reset() {
      revision_ = 0;
      cap_ = init_cap_;
    }
  };

  // two sided bfsが終了した理由
  enum reason_for_finishing_bfs_t {
    kQsIsEmpty,
    kQtIsEmpty,
  };

  bool two_sided_bfs(int s, int t) {
    queue<int> qs, qt;
    qs.push(s); qt.push(t);
    level[s].first = level[t].second = 0;
    bfs_revision[s] = s_side_bfs_revision;
    bfs_revision[t] = t_side_bfs_revision;

    int slevel = 0, tlevel = 0;
    while (sz(qs) != 0 && sz(qt) != 0) {
      if (sz(qs) < sz(qt)) {
        int size = sz(qs);
        FOR(_, size) {
          const int v = qs.front(); qs.pop();
          for (auto t : e[v]) {
            if (t.cap(graph_revision) == 0 || bfs_revision[t.to] == s_side_bfs_revision) continue;
            if (bfs_revision[t.to] == t_side_bfs_revision) {
              return true; // path was found.
            }
            bfs_revision[t.to] = s_side_bfs_revision;
            level[t.to].first = slevel + 1;
            qs.push(t.to);
          }
        }
        slevel++;
      } else {
        int size = sz(qt);
        FOR(_, size) {
          const int v = qt.front(); qt.pop();
          for (auto t : e[v]) {
            if (e[t.to][t.reverse].cap(graph_revision) == 0 || bfs_revision[t.to] == t_side_bfs_revision) continue;
            if (bfs_revision[t.to] == s_side_bfs_revision) {
              return true; // path was found.
            }
            bfs_revision[t.to] = t_side_bfs_revision;
            level[t.to].second = tlevel + 1;
            qt.push(t.to);
          }
        }
      }
      tlevel++;
    }

  reason_for_finishing_bfs = (qs.empty()) ? kQsIsEmpty : kQtIsEmpty;
    return false;
  }

  int dfs(int v, int t, bool use_slevel, int f) {
    if (v == t) return f;
    if (dfs_revision[v] != bfs_revision[v]) {
      dfs_revision[v] = bfs_revision[v];
      iter[v] = 0;
    }
    for (int &i = iter[v]; i < sz(e[v]); i++) {
      E& _e = e[v][i];
      if (_e.cap(graph_revision) == 0 || bfs_revision[_e.to] / 2 != s_side_bfs_revision / 2) continue;

      bool rec;
      if (use_slevel) rec = bfs_revision[_e.to] == t_side_bfs_revision || level[v].first < level[_e.to].first;
      else rec = level[v].second > level[_e.to].second;
      if (!rec) continue;

      bool next_slevel = use_slevel && bfs_revision[_e.to] == s_side_bfs_revision;
      int d = dfs(_e.to, t, next_slevel, min(f, _e.cap(graph_revision)));
      if (d > 0) {
        _e.add_cap(-d, graph_revision);
        e[_e.to][_e.reverse].add_cap(d, graph_revision);
        return d;
      }
    }
    return 0;
  }

  void add_undirected_edge(int f, int t, int c) {
    e[f].push_back(E(t, sz(e[t]), c));
    e[t].push_back(E(f, sz(e[f]) - 1, c));
  }

  void reset_rivision() {
    FOR(v, n) for (auto& e_ : e[v]) e_.reset();
    memset(bfs_revision.data(), 0, sizeof(bfs_revision[0]) * bfs_revision.size());
    memset(dfs_revision.data(), 0, sizeof(dfs_revision[0]) * dfs_revision.size());
    s_side_bfs_revision = 2;
    t_side_bfs_revision = 3;
    graph_revision = 0;
  }

public:
  dinic_twosided() : n(0) {}
  dinic_twosided(const G& g)
    : n(g.num_vertices()), level(n), iter(n), bfs_revision(n), dfs_revision(n), e(n),
    s_side_bfs_revision(2), t_side_bfs_revision(3), graph_revision(0) {
    FOR(v, n) for (auto& e : g.edges(v)) {
      add_undirected_edge(v, to(e), 1);
    }
  }
  dinic_twosided(const vector<pair<V, V>>& edges, int num_vs)
    : n(num_vs), level(n), iter(n), bfs_revision(n), dfs_revision(n), e(n),
    s_side_bfs_revision(2), t_side_bfs_revision(3), graph_revision(0) {
    for (auto& uv : edges) {
      add_undirected_edge(uv.first, uv.second, 1);
    }

  }

  int max_flow(int s, int t) {
    assert(s != t);
    int flow = 0;
    s_side_bfs_revision += 2;
    t_side_bfs_revision += 2;
    for (; ; s_side_bfs_revision += 2, t_side_bfs_revision += 2) {
      bool path_found = two_sided_bfs(s, t);
      if (!path_found) return flow;
      int f;
      while ((f = dfs(s, t, true, numeric_limits<int>::max())) > 0) {
        flow += f;
      }
    }
    return flow;
  }

  bool path_dont_exists_to_t(const int v) const {
    if (reason_for_finishing_bfs == kQsIsEmpty) {
      //sから到達可能な頂点のbfs_revisionには、必ずs_side_bfs_revisionが代入されている
      return bfs_revision[v] == s_side_bfs_revision;
    } else {
      //上の逆
      //tから到達不可能 <=>
      return bfs_revision[v] != t_side_bfs_revision;
    }
  }

  //フローを流す前に実行する
  void reset_graph() {
    graph_revision++;
    if (s_side_bfs_revision >= numeric_limits<decltype(s_side_bfs_revision)>::max() / 2) {
      reset_rivision();
    }
  }

  const int n;
  vector<pair<int, int>> level;
  vector<int> iter;
  vector<int> bfs_revision, dfs_revision;
  vector<vector<E>> e;
  int s_side_bfs_revision, t_side_bfs_revision;
  int graph_revision;
  reason_for_finishing_bfs_t reason_for_finishing_bfs;

};
