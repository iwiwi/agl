#include "easy_cui.h"

DEFINE_int32(num_query, 1000, "");
DEFINE_int64(node_pair_random_seed, 922337203685477583LL, "");

class naive {
  #define sz(c) ((int)c.size())

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
  naive(const G& g) 
  : level(g.num_vertices()), iter(g.num_vertices()), e(g.num_vertices()) {
    for(const auto edge : g.edge_list()){
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

G to_directed_graph(G& g){
  vector<pair<V,V>> ret;
  for(auto& e : g.edge_list()) {
      if(e.first < to(e.second)) ret.emplace_back(e.first,to(e.second));
  }
  return G(ret);
}

class cut_with_degree {
public:
  cut_with_degree(const G& g) : degs(g.num_vertices()) {
    for(V v: make_irange(g.num_vertices())){
      for(auto& e : g.edges(v)) {
        V w = to(e);
        degs[v]++;
        degs[w]++;
      }
    }
  }

  int query(V s,V t) const {
    return min(degs[s],degs[t]);
  }
private:
  vector<int> degs;
};

bool check_min_cut_query(const cut_with_degree& cwd, int s,int t,G& g){
  naive nv(g);
  int naive_w = nv.max_flow(s,t);
  int min_deg = cwd.query(s,t);

  JLOG_ADD_OPEN("query") {
    JLOG_PUT("S",s);
    JLOG_PUT("T",t);
    JLOG_PUT("naive", naive_w);
    JLOG_PUT("min_deg", min_deg);
  }
  if(naive_w != min_deg) {
    fprintf(stderr, "unmatched. (S,T) = (%d,%d), naive = %d, min_deg =%d\n",s,t,naive_w,min_deg);
  }
  return naive_w == min_deg;
}

int main(int argc, char **argv) {
  // JLOG_INIT(&argc, argv); called in "easy_cui_init"
  xorshift64star gen_node(FLAGS_node_pair_random_seed);
  G g = easy_cui_init(argc, argv);
  g = to_directed_graph(g);
  cut_with_degree cwd(g);

  int unmatch = 0;
  for(int counter = 0; counter < FLAGS_num_query; counter++) {
    V s = gen_node() % g.num_vertices();
    V t = gen_node() % (g.num_vertices() - 1);
    if(s <= t) t++;
    if(counter % 100 == 0){
      fprintf(stderr, "count/unmatch/all : %d/%d/%d, \n",counter, unmatch, FLAGS_num_query);
    }
    bool is_matched = check_min_cut_query(cwd, s, t, g);
    if(!is_matched) unmatch++;
  }
  JLOG_PUT("result.all", FLAGS_num_query);
  JLOG_PUT("result.match", (FLAGS_num_query - unmatch));
  JLOG_PUT("result.unmatch", unmatch);
  
  return 0;
}
 