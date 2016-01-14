#pragma once

class dinic_naive {
  struct E {
    int to, rev, cap;
    E(int to, int rev, int cap) : to(to), rev(rev), cap(cap) {}
  };

  void add_undirected_edge(int f, int t, int c) {
    e[f].push_back(E(t, sz(e[t]), c));
    e[t].push_back(E(f, sz(e[f]) - 1, c));
  }

  void bfs(int s) {
    level.assign(level.size(), -1);
    queue<int> q;
    level[s] = 0;
    q.push(s);
    while (!q.empty()) {
      int v = q.front(); q.pop();
      for (auto t : e[v]) {
        if (t.cap > 0 && level[t.to] < 0) {
          level[t.to] = level[v] + 1;
          q.push(t.to);
        }
      }
    }
  }

  int dfs(int v, int t, int f) {
    if (v == t) return f;
    for (int &i = iter[v]; i < sz(e[v]); i++) {
      E& _e = e[v][i];
      if (_e.cap > 0 && level[v] < level[_e.to]) {
        int d = dfs(_e.to, t, min(f, _e.cap));
        if (d > 0) {
          _e.cap -= d;
          e[_e.to][_e.rev].cap += d;
          return d;
        }
      }
    }
    return 0;
  }

  static const int INF = (int)1e8;
public:
  dinic_naive(const G& g) 
  : level(g.num_vertices()), iter(g.num_vertices()), e(g.num_vertices()) {
    for(const auto edge : g.edge_list()){
      add_undirected_edge(edge.first, edge.second, 1);
    }
  }

  dinic_naive(const vector<pair<V, V>>& edges, int num_vs)
  : level(num_vs), iter(num_vs), e(num_vs) {
    for(const auto edge : edges){
      add_undirected_edge(edge.first, edge.second, 1);
    }
  }

  int max_flow(int s, int t) {
    assert(s != t);
    int flow = 0;
    while (true) {
      bfs(s);
      if (level[t] < 0) return flow;
      iter.assign(iter.size(), 0);
      int f;
      while ((f = dfs(s, t, INF)) > 0) {
        flow += f;
      }
    }
  }

  vector<int> level,iter;
  vector<vector<E>> e;
};