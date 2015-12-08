#include <vector>
#include <queue>
#include <cstring>
using namespace std;

#define sz(c) ((int)c.size())

struct E {
  int to, rev, cap;
  E(int to, int rev, int cap) : to(to), rev(rev), cap(cap) {}
};

const int MAX_V = 1000;
int level[MAX_V], iter[MAX_V];
vector<E> e[MAX_V];

void add_undirected_edge(int f, int t, int c) {
  e[f].push_back(E(t, sz(e[t]), c));
  e[t].push_back(E(f, sz(e[f]) - 1, c));
}

void add_directed_edge(int f, int t, int c) {
  e[f].push_back(E(t, sz(e[t]), c));
  e[t].push_back(E(f, sz(e[f]) - 1, 0));
}

void bfs(int s) {
  memset(level, -1, sizeof(level));
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

const int INF = (int)1e8;

int max_flow(int s, int t) {
  int flow = 0;
  while (true) {
    bfs(s);
    if (level[t] < 0) return flow;
    memset(iter, 0, sizeof(iter));
    int f;
    while ((f = dfs(s, t, INF)) > 0) {
      flow += f;
    }
  }
}

#include "easy_cui.h"

void init(G& g){
  V n = g.num_vertices();
  CHECK(n < MAX_V);
  for(V v : make_irange(n)) {
    e[v].clear();
  }
  for(auto edge : g.edge_list()){
    add_undirected_edge(edge.first, edge.second, 1);
  }
}

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);

  V s,t;
  while(cin>> s >> t) {
    init(g);
    int w = max_flow(s,t);
    cout << w << endl;
  }

  return 0;
}
