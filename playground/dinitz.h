#pragma once

namespace agl {
namespace cut_tree_internal {
class dinitz {
  struct E {
    int to_, rev_, cap_;
    E(int to, int rev_, int cap_) : to_(to), rev_(rev_), cap_(cap_) {}
  };

  void bfs(int s) {
    level_.assign(level_.size(), -1);
    queue<int> q;
    level_[s] = 0;
    q.push(s);
    while (!q.empty()) {
      int v = q.front(); q.pop();
      for (auto& t : e_[v]) {
        if (t.cap_ > 0 && level_[t.to_] < 0) {
          level_[t.to_] = level_[v] + 1;
          q.push(t.to_);
        }
      }
    }
  }

  int dfs(int v, int t, int f) {
    if (v == t) return f;
    for (int &i = iter_[v]; i < int(e_[v].size()); i++) {
      E& _e = e_[v][i];
      if (_e.cap_ > 0 && level_[v] < level_[_e.to_]) {
        int d = dfs(_e.to_, t, min(f, _e.cap_));
        if (d > 0) {
          _e.cap_ -= d;
          e_[_e.to_][_e.rev_].cap_ += d;
          return d;
        }
      }
    }
    return 0;
  }

  void add_undirected_edge(int f, int t, int c) {
    e_[f].push_back(E(t, int(e_[t].size()), c));
    e_[t].push_back(E(f, int(e_[f].size()) - 1, c));
  }

  void reset_graph() {
    const int n = g.num_vertices();
    for(int v = 0; v < n; v++) e_[v].clear();
    for(int v = 0; v < n; v++) for (auto& e : g.edges(v)) {
      add_undirected_edge(v, agl::to(e), 1);
    }
  }

public:
  dinitz(const G& g)
    : g(g), level_(g.num_vertices()), iter_(g.num_vertices()), e_(g.num_vertices()) {}

  int max_flow(int s, int t) {
    assert(s != t);
    reset_graph();
    int flow = 0;
    while (true) {
      bfs(s);
      if (level_[t] < 0) return flow;
      iter_.assign(iter_.size(), 0);
      int f;
      while ((f = dfs(s, t, numeric_limits<int>::max())) > 0) {
        flow += f;
      }
    }
  }

  vector<E>& edges(V v) { return e_[v]; }
  V to(const E& e) const { return e.to_; }
  int cap(E& e) { return e.cap_; }
  E& rev(const E& e_in) { return e_[e_in.to_][e_in.rev_]; }

private:
  G g;
  vector<int> level_,iter_;
  vector<vector<E>> e_;
};

} //namespace cut_tree_internal
} //namespace agl