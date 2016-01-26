#pragma once
#include "ma_ordering.h"
#include "dinic_twosided.h"

class dinic_using_ma_ordering {
public:
  enum reason_for_finishing_bfs_t {
    kQsIsEmpty,
    kQtIsEmpty,
  };
private:

  class E {
  private:
    int from_, to_;
    int init_cap_;
    int cap_, revision_;
    int priority_;

  public:
    E(int from, int to, int cap) :
      from_(from), to_(to), init_cap_(cap), cap_(cap), revision_(0), priority_(numeric_limits<int>::max()) {
    }

    const int cap(int from, int currenct_revision) {
      if (revision_ != currenct_revision) {
        revision_ = currenct_revision;
        cap_ = init_cap_;
      }
      getcap_counter++;
      return from == from_ ? cap_ : 2 * init_cap_ - cap_;
    }
    void add_cap(int from, int val, int currenct_revision) {
      if (revision_ != currenct_revision) {
        revision_ = currenct_revision;
        cap_ = init_cap_;
      }
      addcap_counter++;
      if(from == from_) cap_ += val;
      else cap_ -= val;
    }

    int to(int from) const {
      return from == from_ ? to_ : from_;
    }

    int priority() const {
      return priority_;
    }
    void priority(int p) {
      priority_ = p;
    }

    void reset() {
      revision_ = 0;
      cap_ = init_cap_;
    }
  };

  int flow_upper_bound(int s, int t) {
    //各頂点の次数を見て、flowのupper boundを計算する
    return min(sz(eref[s]), sz(eref[t]));
  }

  bool two_sided_bfs(int s, int t, const int st_flow_upper_bound) {
    queue<int> qs, qt;
    qs.push(s); qt.push(t);
    level[s].first = level[t].second = 0;
    bfs_revision[s] = s_side_bfs_revision;
    bfs_revision[t] = t_side_bfs_revision;

    int slevel = 0, tlevel = 0;
    while (sz(qs) != 0 && sz(qt) != 0) {
      bool path_found = false;
      if (sz(qs) < sz(qt)) { //todo : デバッグ用なので直す
        int size = sz(qs);
        FOR(_, size) {
          const int v = qs.front(); qs.pop();
          FOR(itr, sz(eref[v])) {
            auto& t = *eref[v][itr];
            if (t.priority() > st_flow_upper_bound) break;
            const int to = t.to(v);
            if (t.cap(v, graph_revision) == 0 || bfs_revision[to] == s_side_bfs_revision) continue;
            if (bfs_revision[to] == t_side_bfs_revision) {
              path_found = true;
              continue;
            }
            bfs_revision[to] = s_side_bfs_revision;
            level[to].first = slevel + 1;
            qs.push(to);
          }
        }
        slevel++;
      } else {
        int size = sz(qt);
        FOR(_, size) {
          const int v = qt.front(); qt.pop();
          FOR(itr, sz(eref[v])) {
            auto& t = *eref[v][itr];
            if (t.priority() > st_flow_upper_bound) break;
            const int to = t.to(v);
            if (t.cap(to, graph_revision) == 0 || bfs_revision[to] == t_side_bfs_revision) continue;
            if (bfs_revision[to] == s_side_bfs_revision) {
              path_found = true;
              continue;
            }
            bfs_revision[to] = t_side_bfs_revision;
            level[to].second = tlevel + 1;
            qt.push(to);
          }
        }
        tlevel++;
      }
      if (path_found) return true;
    }

    reason_for_finishing_bfs = (qs.empty()) ? kQsIsEmpty : kQtIsEmpty;
    return false;
  }

  int dfs(int v, int t, bool use_slevel, int f,const int st_flow_upper_bound) {
    if (v == t) return f;
    if (dfs_revision[v] != bfs_revision[v]) {
      dfs_revision[v] = bfs_revision[v];
      iter[v] = 0;
    }

    for (int &i = iter[v]; i < sz(eref[v]); i++) {
      E& _e = *eref[v][i];
      if (_e.priority() > st_flow_upper_bound) break;
      const int to = _e.to(v);
      const int cap = _e.cap(v, graph_revision);
      if (cap == 0 || bfs_revision[to] / 2 != s_side_bfs_revision / 2) continue;

      bool rec;
      if (use_slevel) rec = bfs_revision[to] == t_side_bfs_revision || level[v].first < level[to].first;
      else rec = bfs_revision[to] == t_side_bfs_revision && level[v].second > level[to].second;
      if (!rec) continue;

      bool next_slevel = use_slevel && bfs_revision[to] == s_side_bfs_revision;
      int d = dfs(to, t, next_slevel, min(f, cap), st_flow_upper_bound);
      if (d > 0) {
        _e.add_cap(v, -d, graph_revision);
        return d;
      }
    }
    return 0;
  }

  void construct_edge(const vector<pair<V, V>>& edges, int num_vs) {
    for (const auto& edge : edges) {
      pool.emplace_back(edge.first, edge.second, 1);
    }
    // これ以降poolのサイズは変更されないので、E*のアドレスも変わらない -> e[i]にポインタを渡す
    FOR(i,sz(edges)) {
      const auto& edge = edges[i];
      eref[edge.first].push_back(&pool[i]);
      eref[edge.second].push_back(&pool[i]);
    }
  }

  void reset_rivision() {
    FOR(v, n) for (auto& e_ : eref[v]) e_->reset();
    memset(bfs_revision.data(), 0, sizeof(bfs_revision[0]) * bfs_revision.size());
    memset(dfs_revision.data(), 0, sizeof(dfs_revision[0]) * dfs_revision.size());
    s_side_bfs_revision = 2;
    t_side_bfs_revision = 3;
    graph_revision = 0;
  }

public:
  dinic_using_ma_ordering() : n(0) {}

  dinic_using_ma_ordering(const vector<pair<V, V>>& edges, int num_vs)
    : n(num_vs), level(n), iter(n), bfs_revision(n), dfs_revision(n), pool(), eref(n),
    s_side_bfs_revision(2), t_side_bfs_revision(3), graph_revision(0) {

    construct_edge(edges, num_vs);

    vector<V> ma_order = calc_ma_ordering(edges, num_vs);
    vector<int> ma_order_inv(sz(ma_order));
    FOR(i, num_vs) {
      ma_order_inv[ma_order[i]] = i;
    }

    //各辺をma_orderでsortして、priorityを計算する
    FOR(v, num_vs) {
      auto& to_edges = eref[v];
      sort(to_edges.begin(), to_edges.end(), [v, &ma_order_inv](const E* l, const E* r) {
        return ma_order_inv[l->to(v)] < ma_order_inv[r->to(v)];
      });

      // priorityを設定する
      // priority = min( index(edge.from), index(edge.to) );
      FOR(index, sz(to_edges)) {
        E* cur = to_edges[index];
        cur->priority(min(cur->priority(), index));
      }
    }

    //priorityでsortする
    FOR(v, num_vs) {
      auto& to_edges = eref[v];
      sort(to_edges.begin(), to_edges.end(), [](const E* l, const E* r) {
        return l->priority() < r->priority();
      });
    }

    //{
    //  puts("graph G {");
    //  FOR(from, num_vs) {
    //    for (auto& to_edge : eref[from]) {
    //      printf("%d -- %d\n", from, to_edge->to(from));
    //      break;
    //    }
    //  }
    //  puts("}");
    //}
  }

  int max_flow_core(int s, int t) {
    assert(s != t);
    int flow = 0;
    s_side_bfs_revision += 2;
    t_side_bfs_revision += 2;
    const int st_flow_upper_bound = flow_upper_bound(s, t);
    for (; ; s_side_bfs_revision += 2, t_side_bfs_revision += 2) {
      bool path_found = two_sided_bfs(s, t, st_flow_upper_bound);
      if (!path_found)
        return flow;
      while (true) {
        int f = dfs(s, t, true, numeric_limits<int>::max(), st_flow_upper_bound);
        if (f == 0) break;
        flow += f;
      }
    }
    return flow;
  }

  int max_flow(int s, int t) {
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
  vector<E> pool;
  vector<vector<E*>> eref;
  int s_side_bfs_revision, t_side_bfs_revision;
  int graph_revision;
  reason_for_finishing_bfs_t reason_for_finishing_bfs;

};