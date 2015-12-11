#include "easy_cui.h"

DEFINE_int64(solver_iter, 50, "");

namespace naive {
#define sz(c) ((int)c.size())

struct E {
  int to, rev, cap;
  E(int to, int rev, int cap) : to(to), rev(rev), cap(cap) {}
};

const int MAX_V = 5000;
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
  assert(s != t);
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
} // naive

namespace naive {
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
} // naive

class min_cut_query_with_random_contraction {

  int query_dfs(V v, V t, int cost,V par = -1) const {
    if(v == t) return cost;
    for(const auto& to_cost : binary_tree_edges_[v]) {
        V to; int edge_cost; tie(to,edge_cost) = to_cost;
        if(to == par) continue;
        int ncost = query_dfs(to, t, min(cost, edge_cost), v);
        if(ncost != -1) return ncost;
    }

    return -1;
  }

  //O(E uf(E))
  int out_edge_weight(int uf_index) {
    int w = 0;
    for(const auto& edge : initial_edges_) {
      V u = edge.first, v = to(edge.second);
      if(uf_.is_same(uf_index, u) ^ uf_.is_same(uf_index, v)) {
        w++;
      }
    }

    return w;
  }

  void contraction(vector<int>& ancestor, V u, V v) {
      assert(!uf_.is_same(u, v));
      V new_vertex = (V)binary_tree_edges_.size();
      u = uf_.root(u), v = uf_.root(v);
      int uw = out_edge_weight(u), vw = out_edge_weight(v);
      int ua = ancestor[u] , va = ancestor[v];
      uf_.unite(u, v);
      ((uf_.root(u) == u) ? ancestor[u] : ancestor[v]) = new_vertex;

      // 縮約した結果を元にbinary treeの構築を進める
      binary_tree_edges_.emplace_back(); // new_vertex分の確保
      binary_tree_edges_[new_vertex].emplace_back(ua, uw);
      binary_tree_edges_[ua].emplace_back(new_vertex, uw);
      binary_tree_edges_[new_vertex].emplace_back(va, vw);
      binary_tree_edges_[va].emplace_back(new_vertex, vw);       
    };

 public:
  // g の辺を使い、グラフを作成する(勝手にundirectedとして読み替えている)
  min_cut_query_with_random_contraction(G& g) :
    num_vertices_(g.num_vertices()),
    initial_edges_(g.edge_list()),
    uf_(g.num_vertices()),
    binary_tree_edges_(g.num_vertices()) {
    //unordered -> weight = 1なので、random_shuffleでよい
    //重みがあるならBITで。(stochastic acceptanceはupdateありだと使えない)
    std::random_shuffle(initial_edges_.begin(),initial_edges_.end());
    vector<int> ancestor(g.num_vertices());
    std::iota(ancestor.begin(),ancestor.end(), 0);

    // 次数1の頂点を先に縮約する
    {
      vector<int> degree;
      queue<int> q;
      for(V v : make_irange(g.num_vertices())) {
        degree.push_back(g.degree(v, kFwd) + g.degree(v, kBwd));
        if(degree[v] == 1) q.push(v);
      }
      while(!q.empty()) {
        V v = q.front(); q.pop();
        for(D dir : direction_range()) {
          for(auto edge : g.neighbors(v, dir)) {
            V u = to(edge);
            degree[u]--;
            if(degree[u] >= 1) contraction(ancestor, v, u); //つながってる先へ縮約
            if(degree[u] == 1) q.push(u);
          }
        }
      }
    }

    //O(VE uf(E))
    for(const auto& edge : initial_edges_) {
      V u = edge.first, v = to(edge.second);
      if(uf_.is_same(u,v)) continue;
      //縮約
      contraction(ancestor, u, v);
    }
    // gが連結とは限らないので、uv_costs.size() == g.num_vertices() - 1ではない
    assert(int(binary_tree_edges_.size() - g.num_vertices()) <= g.num_vertices() - 1);
  }

  void debug_output_graph(){
    ostringstream oss;
    for(V v = num_vertices_; v < (V)binary_tree_edges_.size(); v++) {
      for(auto& to_cost : binary_tree_edges_[v]) {
        V to; int edge_cost; tie(to,edge_cost) = to_cost;
        if(v < to) continue;
        oss << v << " " << to << endl;
      }
    }
    istringstream iss(oss.str());
    auto g = read_graph_tsv<G>(iss);
    graphviz gv(g);
    for(V v = num_vertices_; v < (V)binary_tree_edges_.size(); v++) {
      for(auto& to_cost : binary_tree_edges_[v]) {
        V to; int edge_cost; tie(to,edge_cost) = to_cost;
        if(v < to) continue;
        gv.set_edge_property(v, to, "weight", to_string(edge_cost));
      }
    }
    gv.generate_png("binary.png");
  }

  int query(V u,V v) const {
    CHECK(u != v);
    CHECK(u < num_vertices_ && v < num_vertices_);
    int ans = query_dfs(u,v, numeric_limits<int>::max());
    if(ans == -1) {
      return 0; // 到達できなかった
    }
    return ans;
  }

private:
  const int num_vertices_;
  typename G::edge_list_type initial_edges_;
  union_find uf_;
  vector<vector<pair<V,int>>> binary_tree_edges_;
};

class min_cut_query {
 public:
  min_cut_query(G& g) {
    for(int i = 0; i < FLAGS_solver_iter; i++) {
      solvers_.emplace_back(g);
    }
  }

  int query(int u,int v) const {
    int ans = numeric_limits<int>::max();
    for(const auto& solver : solvers_){
      ans = min(ans, solver.query(u, v));
    }
    return ans;
  }
 private:
  vector<min_cut_query_with_random_contraction> solvers_;
};


G to_directed_graph(G& g){
  vector<pair<V,V>> ret;
  for(auto& e : g.edge_list()) {
      if(e.first < to(e.second)) ret.emplace_back(e.first,to(e.second));
  }
  return G(ret);
}

bool check_min_cut_query(const min_cut_query& mcq,int s,int t,G& g){
  naive::init(g);
  int naive_w = naive::max_flow(s,t);
  int mcq_w = mcq.query(s,t);
  if(naive_w == mcq_w) return true; //matched
 
  CHECK(naive_w <= mcq_w);
  JLOG_ADD_OPEN("unmatched") {
    JLOG_PUT("S",s);
    JLOG_PUT("T",t);
    JLOG_PUT("naive", naive_w);
    JLOG_PUT("mcq_w", mcq_w);
  }
  return false;
}

int main(int argc, char **argv) {
  // JLOG_INIT(&argc, argv); called in "easy_cui_init"
  G g = easy_cui_init(argc, argv);
  g = to_directed_graph(g);
  min_cut_query mcq(g);
  // mcq.debug_output_graph();

  JLOG_PUT("argv.mincut_tree_iter", FLAGS_solver_iter);

  const int all = g.num_vertices() * (g.num_vertices() - 1) / 2;
  int counter = 0,unmatch = 0;
  JLOG_ADD_OPEN("sampleing") {
      for(V s = 0; s < g.num_vertices(); s++) {
        for(V t = s + 1; t < g.num_vertices(); t++) {
          if(counter % 100 == 0){
            fprintf(stderr, "count/unmatch/all : %d/%d/%d, \n",counter, unmatch, all);
          }
          bool is_matched = check_min_cut_query(mcq, s, t, g);
          if(!is_matched) unmatch++;
          counter++;
        }
      }
  }
  JLOG_PUT("result.all", all);
  JLOG_PUT("result.match", (all - unmatch));
  JLOG_PUT("result.unmatch", unmatch);

  return 0;
}
 