#include "easy_cui.h"
DEFINE_int32(num_query, 1000, "");
DEFINE_int64(node_pair_random_seed, 922337203685477583LL, "");
DEFINE_string(exec,"both", "dinic, push_relabel, both");

#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())

class dinic {

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
  dinic(const G& g)
   : level(g.num_vertices()), iter(g.num_vertices()), e(g.num_vertices()) {
   for (const auto edge : g.edge_list()) {
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

  vector<int> level, iter;
  vector<vector<E>> e;
};

class push_relabel {
  struct E {
    int to, rev, cap;
    E(int to, int rev, int cap) : to(to), rev(rev), cap(cap) {}
  };

  void initialize_preflow(int s) {
    h[s] = n;
    for (auto& fwd : edges[s]) {
      auto& bwd = edges[fwd.to][fwd.rev];
      int preflow = fwd.cap;
      excess[s] -= preflow;
      excess[fwd.to] += preflow;
      fwd.cap -= preflow;
      bwd.cap += preflow;
    }
  }

  void push(int u, E& fwd) {
    auto& bwd = edges[fwd.to][fwd.rev];
    int preflow = min(excess[u], fwd.cap);
    excess[u] -= preflow;
    excess[fwd.to] += preflow;
    fwd.cap -= preflow;
    bwd.cap += preflow;
  }

public:
  void add_undirected_edge(int f, int t, int c) {
    edges[f].push_back(E(t, sz(edges[t]), c));
    edges[t].push_back(E(f, sz(edges[f]) - 1, c));
  }

  push_relabel(const G& g) : n(g.num_vertices()), edges(g.num_vertices()), excess(g.num_vertices()), h(g.num_vertices()) {
    for (const auto edge : g.edge_list()) {
     add_undirected_edge(edge.first, edge.second, 1);
    }
  }

  int max_flow(int s, int t) {
    queue<int> q;
    vector<char> in_queue(n);
    in_queue[s] = in_queue[t] = true; //s,tはqueueに入らないようにする
  
    initialize_preflow(s);
    for (auto& e : edges[s]) {
      if(e.to == t) continue;
      if (excess[e.to] > 0) {
        q.push(e.to);
        in_queue[e.to] = true;
      }
    }

    while (!q.empty()) {
      // discharge
      int u = q.front();
      assert(excess[u] > 0);
      int min_h = 2 * n; //"h" is always less than 2 * n
      for (auto& e : edges[u]) {
        if (e.cap == 0) continue;
        if (h[u] == h[e.to] + 1) {
          push(u, e);
          if (!in_queue[e.to]) {
            in_queue[e.to] = true;
            q.push(e.to);
          } 
          if(excess[u] == 0) break; // 継続しない
        } else { // h[u] <= h[e.to]
          min_h = min(min_h, h[e.to]);
        }
      }
      if (excess[u] == 0) {
        in_queue[u] = false;
        q.pop();
      } else { 
        // relabel
        assert(h[u] < min_h + 1);
        h[u] = min_h + 1; // 流し足りていない -> 高くする
        // queueからpopしないので、whileに戻りもう一度同じ頂点から操作が開始される
      }
    }

    return excess[t];
  }

  const int n;
  vector<vector<E>> edges;
  vector<int> excess;
  vector<int> h;
};

G to_directed_graph(G& g){
  vector<pair<V,V>> ret;
  for(auto& e : g.edge_list()) {
      if(e.first < to(e.second)) ret.emplace_back(e.first,to(e.second));
  }
  return G(ret);
}

bool check_min_cut_query(V s,V t,G& g) {
  if(FLAGS_exec == "dinic") {
    dinic dc(g);
    int dc_flow = dc.max_flow(s, t);
    return true;
  } else if(FLAGS_exec == "push_relabel") {
    push_relabel pr(g);
    int pr_flow = pr.max_flow(s, t);
    return true;
  } else {
    dinic dc(g);
    push_relabel pr(g);
    int dc_flow = dc.max_flow(s, t);
    int pr_flow = pr.max_flow(s, t);
    return dc_flow == pr_flow;
  }
}


int main(int argc, char** argv) {
  xorshift64star gen_node(FLAGS_node_pair_random_seed);
  G g = easy_cui_init(argc, argv);
  g = to_directed_graph(g);
  push_relabel pr(g);

  int unmatch = 0;
  for(int counter = 0; counter < FLAGS_num_query; counter++) {
    V s = gen_node() % g.num_vertices();
    V t = gen_node() % (g.num_vertices() - 1);
    if(s <= t) t++;
    if(counter % 1 == 0){
      fprintf(stderr, "count/unmatch/all : %d/%d/%d, \n",counter, unmatch, FLAGS_num_query);
    }
    bool is_matched = check_min_cut_query(s, t, g);
    if(!is_matched) unmatch++;
  }
  JLOG_PUT("result.all", FLAGS_num_query);
  JLOG_PUT("result.match", (FLAGS_num_query - unmatch));
  JLOG_PUT("result.unmatch", unmatch);


}