#include "easy_cui.h"

DEFINE_bool(prune_if_degree_eq_1, true, "");
DEFINE_int32(solver_iter, 50, "");
DEFINE_int32(num_query, 1000, "");
DEFINE_int64(node_pair_random_seed, 922337203685477583LL, "");

#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())

#include "dinic_twosided.h"

class naive {

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

  void build_depth() {
    depth_.resize(sz(binary_tree_edges_), -1);
    parent_weight_.resize(sz(binary_tree_edges_), make_pair(-2, -2));

    FOR(v, num_vertices_) {
      if(depth_[v] >= 0) continue;

      depth_[v] = 0;
      parent_weight_[v] = make_pair(-1, -1);
      queue<int> q;
      q.push(v);
      while (!q.empty()) {
        int u = q.front(); q.pop();
        for(auto& to_weight : binary_tree_edges_[u]) {
          int to, weight; tie(to,weight) = to_weight;
          if(depth_[to] >= 0) continue;
          depth_[to] = depth_[u] + 1;
          parent_weight_[to].first = u;
          parent_weight_[to].second = weight;
          // printf("%d - %d\n",u, to);
          q.push(to);
        }
      }
    }
  }

  void contraction(vector<int>& ancestor,vector<unordered_set<pair<V,V>>>& contraction_graph_edges , V u, V v) {
      assert(!uf_.is_same(u, v));
      V new_vertex = (V)binary_tree_edges_.size();
      u = uf_.root(u), v = uf_.root(v);
      int uw = (int)contraction_graph_edges[u].size(), vw = (int)contraction_graph_edges[v].size();
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
      auto& uset = contraction_graph_edges[u];
      auto& vset = contraction_graph_edges[v];
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
    uf_(g.num_vertices()),
    binary_tree_edges_(g.num_vertices()) {
    //unordered -> weight = 1なので、shuffleでよい
    //重みがあるならBITで。(stochastic acceptanceはupdateありだと使えない)
    vector<unordered_set<pair<V,V>>> contraction_graph_edges(g.num_vertices());
    typename G::edge_list_type initial_edges(g.edge_list());
    std::shuffle(initial_edges.begin(),initial_edges.end(), agl::random);
    vector<int> ancestor(g.num_vertices());
    std::iota(ancestor.begin(),ancestor.end(), 0);

    for(auto edge : initial_edges) {
      V u = edge.first, v = to(edge.second);
      contraction_graph_edges[u].emplace(u, v);
      contraction_graph_edges[v].emplace(v, u);
    }

    // 次数1の頂点を先に縮約する
    if(FLAGS_prune_if_degree_eq_1){
      cerr << "First, prune node which degree eq 1." << endl;
      int pruned = 0;
      queue<int> q;
      for(V v : make_irange(g.num_vertices())) {
        if(contraction_graph_edges[v].size() == 1) q.push(v);
      }
      while(!q.empty()) {
        V v = q.front(); q.pop();
        pruned++;
        V u;
        if(contraction_graph_edges[v].size() == 0) continue;
        tie(std::ignore, u) = *contraction_graph_edges[v].begin();
        contraction(ancestor, contraction_graph_edges, u, v);
        if(contraction_graph_edges[u].size() == 1) q.push(u);
      }
      cerr << pruned << " node(s) was pruned." << endl;
    }

    //O(E log(E))
    for(const auto& edge : initial_edges) {
      V u = edge.first, v = to(edge.second);
      if(uf_.is_same(u,v)) continue;
        contraction(ancestor, contraction_graph_edges, u, v);
    }
    // gが連結とは限らないので、uv_costs.size() == g.num_vertices() - 1ではない
    assert(int(binary_tree_edges_.size() - g.num_vertices()) <= g.num_vertices() - 1);

    build_depth();
  }

    int query(V u, V v) const {
    CHECK(u != v);
    CHECK(u < num_vertices_ && v < num_vertices_);
    int ans = numeric_limits<int>::max();
    while (u != v) {
      if(u == -1 || v == -1) return 0; // disconnect
      if (depth_[u] > depth_[v]) {
        ans = min(ans, parent_weight_[u].second);
        u = parent_weight_[u].first;
      } else {
        ans = min(ans, parent_weight_[v].second);
        v = parent_weight_[v].first;
      }
    }
    return ans;
  }

  void single_source_mincut_dfs(V v,V par,int mn, vector<int>& ans) const {
    if(v < num_vertices_) {
      ans[v] = mn;
    }
    for(auto& to : binary_tree_edges_[v]) {
      if(to.first == par) continue;
      int nmn = min(mn, to.second);
      single_source_mincut_dfs(to.first, v, nmn, ans);
    }
  }


  vector<int> single_source_mincut(V u) const {
    vector<int> ans(num_vertices_, numeric_limits<int>::max());
    single_source_mincut_dfs(u,-1,numeric_limits<int>::max(), ans);
    return ans;
  }

private:
  const int num_vertices_;
  union_find uf_;

  vector<vector<pair<V,int>>> binary_tree_edges_;
  vector<pair<V, int>> parent_weight_;
  vector<int> depth_;
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

  vector<int> single_source_mincut(V v) const {
    vector<int> ret;
    for(const auto& solver : solvers_){
      auto cur_ans = solver.single_source_mincut(v);
      if(sz(ret) == 0) ret = move(cur_ans);
      else {
        FOR(i,sz(ret)) ret[i] = min(ret[i], cur_ans[i]);
      }
    }
    return ret;
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

DEFINE_string(method,"test","test, single_source_mincut");

void test(G&& g) {
    min_cut_query mcq(g);
  // mcq.debug_output_graph();

  JLOG_PUT("argv.mincut_tree_iter", FLAGS_solver_iter);

  xorshift64star gen_node(FLAGS_node_pair_random_seed);
  int unmatch = 0;
  for(int counter = 0; counter < FLAGS_num_query; counter++) {
    V s = gen_node() % g.num_vertices();
    V t = gen_node() % (g.num_vertices() - 1);
    if(s <= t) t++;
    if(counter % 100 == 0){
      fprintf(stderr, "count/unmatch/all : %d/%d/%d, \n",counter, unmatch, FLAGS_num_query);
    }
    bool is_matched = check_min_cut_query(mcq, s, t, g);
    if(!is_matched) unmatch++;
  }
  JLOG_PUT("result.all", FLAGS_num_query);
  JLOG_PUT("result.match", (FLAGS_num_query - unmatch));
  JLOG_PUT("result.unmatch", unmatch);
}

string graph_name() {
  string x = FLAGS_graph;
  string ret;
  for(int i = sz(x) - 1; i >= 0; i--) {
   if(x[i] == '/' || x[i] == '\\') break;
    ret.push_back(x[i]);
  }
  reverse(ret.begin(),ret.end());
  return ret;
}

DEFINE_string(single_source_mincut_output, "random_contraction_ssm.data", "");
void single_source_mincut(G&& g) {
    min_cut_query mcq(g);
    size_t max_deg = 0;
    int max_deg_v = -1;
    FOR(v, g.num_vertices()) {
      if(g.degree(v, D(0)) + g.degree(D(1)) > max_deg) {
        max_deg = g.degree(v, D(0)) + g.degree(D(1));
        max_deg_v = v;
      }
    }

    auto anses = mcq.single_source_mincut(max_deg_v);
    FILE* fp = fopen(FLAGS_single_source_mincut_output.c_str(), "w");
    FOR(v,g.num_vertices()) {
      if(v == max_deg_v) continue;
      
      fprintf(fp, "%d %d %d\n",max_deg_v, v, anses[v]);
    }  
    fclose(fp);
}

void single_source_mincut_with_deg(G&& g) {
    size_t max_deg = 0;
    int max_deg_v = -1;
    FOR(v, g.num_vertices()) {
      if(g.degree(v) > max_deg) {
        max_deg = g.degree(v);
        max_deg_v = v;
      }
    }

    FILE* fp = fopen(FLAGS_single_source_mincut_output.c_str(), "w");
    FOR(v,g.num_vertices()) {
      if(v == max_deg_v) continue;
      int xx = g.degree(v, D(0)) + g.degree(v, D(1));
      fprintf(fp, "%d %d %d\n",max_deg_v, v, xx);
    }  
    fclose(fp);
}


int main(int argc, char **argv) {
  // JLOG_INIT(&argc, argv); called in "easy_cui_init"
  G g = easy_cui_init(argc, argv);
  g = to_directed_graph(g);

  if(FLAGS_method == "test") {
    test(move(g));
  } else if(FLAGS_method == "single_source_mincut") {
    single_source_mincut(move(g));
  } else if(FLAGS_method == "single_source_mincut_with_deg") {
    single_source_mincut_with_deg(move(g));
  } else {
      fprintf(stderr, "unrecognized option '-method=%s'\n", FLAGS_method.c_str());
      exit(-1);
  }

  
  return 0;
}
 