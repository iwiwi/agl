#pragma once

DEFINE_int32(flow_iter, 1, "");

long long getcap_counter = 0;
long long addcap_counter = 0;

class dinic_twosided {
public:
  // two sided bfsが終了した理由
  enum reason_for_finishing_bfs_t {
    kQsIsEmpty,
    kQtIsEmpty,
  };
private:

  class E {
    int init_cap_;
    int cap_, revision_;
  public:
    E(int to, int reverse, int cap) :
      init_cap_(cap), cap_(cap), revision_(0), to(to), reverse(reverse) {
    }

    int to, reverse;
    const int cap(int currenct_revision) {
      if (revision_ != currenct_revision) {
        revision_ = currenct_revision;
        cap_ = init_cap_;
      }
      getcap_counter++;
      return cap_;
    }
    const int cap_without_revision() const {
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

  class vecE {
  public:
    vecE() : rem_size_(0) {}

    void add(E&& to) {
      to_.push_back(to);
      rem_size_++;
    }

    bool match_revision(int revision) {
      return revision_ == revision;
    }

    void reset_graph(int revision) {
      rem_size_ = sz(to_);
      revision_ = revision;
    }

    void reset_revision() {
      rem_size_ = sz(to_);
      for (auto& e : to_) e.reset();
      revision_ = 0;
    }

    bool empty() const { return rem_size_ == 0; }
    int size() const { return rem_size_; }


    E& operator[](const int id) { return to_[id]; };

    int rem_size_, revision_;
    vector<E> to_;
  };

  bool two_sided_bfs(int s, int t) {
    queue<int> qs, qt;
    qs.push(s); qt.push(t);
    level[s].first = level[t].second = 0;
    bfs_revision[s] = s_side_bfs_revision;
    bfs_revision[t] = t_side_bfs_revision;

    int slevel = 0, tlevel = 0;
    while (sz(qs) != 0 && sz(qt) != 0) {
      bool path_found = false;
      if (sz(qs) < sz(qt)) {
        int size = sz(qs);
        FOR(_, size) {
          const int v = qs.front(); qs.pop();
          if (!e[v].match_revision(graph_revision)) {
            e[v].reset_graph(graph_revision);
          }
          const int esize = e[v].size();
          FOR(i, esize) {
            auto& t = e[v][i];
            if (bfs_revision[t.to] == s_side_bfs_revision) continue;
            if (bfs_revision[t.to] == t_side_bfs_revision) {
              path_found = true;
              continue;
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
          if (!e[v].match_revision(graph_revision)) {
            e[v].reset_graph(graph_revision);
          }
          const int esize = sz(e[v].to_);
          FOR(i, esize) {
            auto& t = e[v][i];
            if (bfs_revision[t.to] == t_side_bfs_revision || e[t.to][t.reverse].cap(graph_revision) == 0) continue;
            if (bfs_revision[t.to] == s_side_bfs_revision) {
              path_found = true;
              continue;
            }
            bfs_revision[t.to] = t_side_bfs_revision;
            level[t.to].second = tlevel + 1;
            qt.push(t.to);
          }
        }
        tlevel++;
      }
      if (path_found) return true;
    }

    reason_for_finishing_bfs = (qs.empty()) ? kQsIsEmpty : kQtIsEmpty;
    return false;
  }

  void debug() {
    FOR(v, n) {
      if (!e[v].match_revision(graph_revision)) {
        e[v].reset_graph(graph_revision);
      }
      FOR(i,sz(e[v].to_) ) {
        auto& e_ = e[v][i];
        auto& rev = e[e_.to][e_.reverse];
        auto& e2 = e[rev.to][rev.reverse];
        CHECK(&e_ == &e2);

        if (i < e[v].rem_size_) {
          CHECK(e_.cap(graph_revision) > 0);
        } else {
          CHECK(e_.cap(graph_revision) == 0);
        }
      }
    }
  }

  //辺のcapが空っぽになる
  void to_empty(int from, int edge_id) {
    // 逆辺の整合性を保ちつつ辺を入れ替える
    E& cur = e[from][edge_id];
    int swap_id = e[from].rem_size_ - 1;
    if (swap_id != edge_id) {
      E& swaped = e[from][swap_id];
      CHECK(edge_id < e[from].rem_size_ && cur.cap_without_revision() == 0); // cap が空だけど、cap==0の領域に入っていない
      CHECK(swaped.cap(graph_revision) > 0);

      E& cur_rev = e[cur.to][cur.reverse];
      E& swap_rev = e[swaped.to][swaped.reverse];

      std::swap(cur_rev.reverse, swap_rev.reverse);
      std::swap(cur, swaped);
    }

    e[from].rem_size_--;

    // debug();
  }

  //辺のcapが空じゃなくなる
  void to_not_empty(int from, int edge_id) {
    E& cur = e[from][edge_id];
    int swap_id = e[from].rem_size_;
    if (swap_id != edge_id) {
      E& swaped = e[from][swap_id];
      CHECK(edge_id >= e[from].rem_size_ && cur.cap_without_revision() > 0); // cap が空じゃないのに、cap==0の領域に入っている
      CHECK(swaped.cap(graph_revision) == 0);

      E& cur_rev = e[cur.to][cur.reverse];
      E& swap_rev = e[swaped.to][swaped.reverse];

      std::swap(cur_rev.reverse, swap_rev.reverse);
      std::swap(cur, swaped);
    }

    e[from].rem_size_++;

    // debug();
  }


  int dfs(int v, int t, bool use_slevel, int f) {
    if (v == t) return f;
    if (dfs_revision[v] != bfs_revision[v]) {
      dfs_revision[v] = bfs_revision[v];
      iter[v] = 0;
    }
    if (!e[v].match_revision(graph_revision)) {
      e[v].reset_graph(graph_revision);
    }
    for (int &i = iter[v]; i < e[v].size(); i++) {
      E& _e = e[v][i];
      if (bfs_revision[_e.to] / 2 != s_side_bfs_revision / 2) continue;
      const int cap = _e.cap(graph_revision);
      CHECK(cap > 0);

      bool rec;
      if (use_slevel) rec = bfs_revision[_e.to] == t_side_bfs_revision || level[v].first < level[_e.to].first;
      else rec = bfs_revision[_e.to] == t_side_bfs_revision && level[v].second > level[_e.to].second;
      if (!rec) continue;

      bool next_slevel = use_slevel && bfs_revision[_e.to] == s_side_bfs_revision;
      int d = dfs(_e.to, t, next_slevel, min(f, cap));
      if (d > 0) {
        _e.add_cap(-d, graph_revision);
        auto& rev = e[_e.to][_e.reverse];
        rev.add_cap(d, graph_revision);

        if (_e.cap(graph_revision) == 0) {
          //_e.cap が d -> 0 に変化した
          to_empty(v, i);
        }
        if (rev.cap(graph_revision) == d) {
          //rev.cap が 0 -> d に変化した
          to_not_empty(_e.to, _e.reverse);
        }

        return d;
      }
    }
    return 0;
  }

  void add_undirected_edge(int f, int t, int c) {
    e[f].add(E(t, sz(e[t]), c));
    e[t].add(E(f, sz(e[f]) - 1, c));
  }

  void reset_rivision() {
    FOR(v, n) e[v].reset_revision();
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

  int max_flow_core(int s, int t) {
    assert(s != t);
    int flow = 0;
    s_side_bfs_revision += 2;
    t_side_bfs_revision += 2;
    for (; ; s_side_bfs_revision += 2, t_side_bfs_revision += 2) {
      bool path_found = two_sided_bfs(s, t);
      if (!path_found) return flow;
      while (true) {
        int f = dfs(s, t, true, numeric_limits<int>::max());
        if (f == 0) break;
        flow += f;
      }
    }
    return flow;
  }

  int max_flow(int s, int t) {
    // cout << s << " " << t << endl;
    auto b1 = getcap_counter;
    int ans = max_flow_core(s, t);
    auto b2 = getcap_counter;
    FOR(_, FLAGS_flow_iter - 1) {
      reset_graph();
      max_flow_core(s, t);
      auto b3 = getcap_counter;
      CHECK(b3 - b2 == b2 - b1);
      b1 = b2;
      b2 = b3;
    }
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

  const int n;
  vector<pair<int, int>> level;
  vector<int> iter;
  vector<int> bfs_revision, dfs_revision;
  vector<vecE> e;
  int s_side_bfs_revision, t_side_bfs_revision;
  int graph_revision;
  reason_for_finishing_bfs_t reason_for_finishing_bfs;

};
