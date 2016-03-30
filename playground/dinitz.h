#pragma once

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

  static const int INF = (int)1e8;
public:
  dinitz(int num_vs)
  : level_(num_vs), iter_(num_vs), e_(num_vs) {
  }

  void add_undirected_edge(int f, int t, int c) {
    e_[f].push_back(E(t, int(e_[t].size()), c));
    e_[t].push_back(E(f, int(e_[f].size()) - 1, c));
  }

  int max_flow(int s, int t) {
    assert(s != t);
    int flow = 0;
    while (true) {
      bfs(s);
      if (level_[t] < 0) return flow;
      iter_.assign(iter_.size(), 0);
      int f;
      while ((f = dfs(s, t, INF)) > 0) {
        flow += f;
      }
    }
  }

  vector<E>& edges(V v) { return e_[v]; }

private:
  vector<int> level_,iter_;
  vector<vector<E>> e_;
};