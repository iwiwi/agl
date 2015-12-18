#include "easy_cui.h"

DEFINE_int32(solver_iter, 50, "");
DEFINE_int32(num_query, 1000, "");

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

  void contraction(vector<int>& ancestor, V u, V v) {
      assert(!uf_.is_same(u, v));
      V new_vertex = (V)binary_tree_edges_.size();
      u = uf_.root(u), v = uf_.root(v);
      int uw = (int)contraction_graph_edges_[u].size(), vw = (int)contraction_graph_edges_[v].size();
      int ua = ancestor[u] , va = ancestor[v];
      uf_.unite(u, v);
      if(uf_.root(u) != u) swap(u, v);
      ancestor[u] = new_vertex;

      // 縮約した結果を元にbinary treeの構築を進める
      binary_tree_edges_.emplace_back(); // new_vertex分の確保
      binary_tree_edges_[new_vertex].emplace_back(ua, uw);
      binary_tree_edges_[ua].emplace_back(new_vertex, uw);
      binary_tree_edges_[new_vertex].emplace_back(va, vw);
      binary_tree_edges_[va].emplace_back(new_vertex, vw);

      //縮約した頂点の辺をまとめる
      auto& uset = contraction_graph_edges_[u];
      auto& vset = contraction_graph_edges_[v];
      if(uset.size() < vset.size()) uset.swap(vset);
      for(const auto& edge : vset) {
        V from, to; tie(from, to) = edge;
        if(uf_.is_same(from, to)){
          auto it = uset.find(make_pair(to, from));
          assert(it != uset.end());
          uset.erase(it);
        } else {
          uset.insert(edge);
        }
      }
      vset.clear();
    }

 public:
  // g の辺を使い、グラフを作成する(勝手にundirectedとして読み替えている)
  min_cut_query_with_random_contraction(G& g) :
    num_vertices_(g.num_vertices()),
    contraction_graph_edges_(g.num_vertices()),
    uf_(g.num_vertices()),
    binary_tree_edges_(g.num_vertices()) {
    //unordered -> weight = 1なので、shuffleでよい
    //重みがあるならBITで。(stochastic acceptanceはupdateありだと使えない)
    typename G::edge_list_type initial_edges(g.edge_list());
    std::shuffle(initial_edges.begin(),initial_edges.end(), agl::random);
    vector<int> ancestor(g.num_vertices());
    std::iota(ancestor.begin(),ancestor.end(), 0);

    for(auto edge : initial_edges) {
      V u = edge.first, v = to(edge.second);
      contraction_graph_edges_[u].emplace(u, v);
      contraction_graph_edges_[v].emplace(v, u);
    }

    // 次数1の頂点を先に縮約する
    {
      queue<int> q;
      for(V v : make_irange(g.num_vertices())) {
        if(contraction_graph_edges_[v].size() == 1) q.push(v);
      }
      while(!q.empty()) {
        V v = q.front(); q.pop();
        V u;
        if(contraction_graph_edges_[v].size() == 0) continue;
        tie(std::ignore, u) = *contraction_graph_edges_[v].begin();
        contraction(ancestor, v, u);
        if(contraction_graph_edges_[u].size() == 1) q.push(u);
      }
    }

    //O(E log(E))
    for(const auto& edge : initial_edges) {
      V u = edge.first, v = to(edge.second);
      if(uf_.is_same(u,v)) continue;
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
  vector<unordered_set<pair<V,V>>> contraction_graph_edges_;
  union_find uf_;
  vector<vector<pair<V,int>>> binary_tree_edges_;
};

class min_cut_query {
 public:
  min_cut_query(G& g) {
    fprintf(stderr, "initialize mincut-query tree ...\n");
    for(int i = 0; i < FLAGS_solver_iter; i++) {
      if(i % 10 == 0) fprintf(stderr, "trees .. %d/%d\n",i,FLAGS_solver_iter);
      solvers_.emplace_back(g);
    }
    fprintf(stderr, "completed.\n");
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
  naive nv(g);
  int naive_w = nv.max_flow(s,t);
  int mcq_w = mcq.query(s,t);
  CHECK(naive_w <= mcq_w);

  JLOG_ADD_OPEN("query") {
    JLOG_PUT("S",s);
    JLOG_PUT("T",t);
    JLOG_PUT("naive", naive_w);
    JLOG_PUT("mcq_w", mcq_w);
  }
  if(naive_w != mcq_w) {
    fprintf(stderr, "unmatched. (S,T) = (%d,%d), naive = %d, mcq_w =%d\n",s,t,naive_w,mcq_w);
  }
  return naive_w == mcq_w;
}

int main(int argc, char **argv) {
  // JLOG_INIT(&argc, argv); called in "easy_cui_init"
  G g = easy_cui_init(argc, argv);
  g = to_directed_graph(g);
  min_cut_query mcq(g);
  // mcq.debug_output_graph();

  JLOG_PUT("argv.mincut_tree_iter", FLAGS_solver_iter);

  int unmatch = 0;
  JLOG_ADD_OPEN("sampleing") {
    for(int counter = 0; counter < FLAGS_num_query; counter++) {
      V s = agl::random() % g.num_vertices();
      V t = agl::random() % (g.num_vertices() - 1);
      if(s <= t) t++;
      if(counter % 100 == 0){
        fprintf(stderr, "count/unmatch/all : %d/%d/%d, \n",counter, unmatch, FLAGS_num_query);
      }
      bool is_matched = check_min_cut_query(mcq, s, t, g);
      if(!is_matched) unmatch++;
      counter++;
    }
  }
  JLOG_PUT("result.all", FLAGS_num_query);
  JLOG_PUT("result.match", (FLAGS_num_query - unmatch));
  JLOG_PUT("result.unmatch", unmatch);
  
  return 0;
}
 