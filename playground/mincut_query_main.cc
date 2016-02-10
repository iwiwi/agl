#include "easy_cui.h"

#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())

#include "dinic_naive.h"
#include "dinic_twosided.h"
#include "random_contraction.h"

DEFINE_int32(num_query, 1000, "");
DEFINE_int64(node_pair_random_seed, 922337203685477583LL, "");
DEFINE_string(method, "test", "test, single_source_mincut");

bool check_min_cut_query(const min_cut_query& mcq, dinic_twosided& dc, int s, int t, G& g) {
  dc.reset_graph();
  int naive_w = dc.max_flow(s, t);
  int mcq_w = mcq.query(s, t);
  if (naive_w > mcq_w) {
    dinic_naive dn(g);
    int val = dn.max_flow(s, t);
    printf("naive = %d, mcq = %d, dinic_twosided = %d\n",val,mcq_w, naive_w);
  }
  CHECK(naive_w <= mcq_w);

  JLOG_ADD_OPEN("query") {
    JLOG_PUT("S", s);
    JLOG_PUT("T", t);
    JLOG_PUT("naive", naive_w);
    JLOG_PUT("mcq_w", mcq_w);
  }
  if (naive_w != mcq_w) {
    fprintf(stderr, "unmatched. (S,T) = (%d,%d), naive = %d, mcq_w =%d\n", s, t, naive_w, mcq_w);
  }
  return naive_w == mcq_w;
}

class random_contraction2 {
  struct degree_weight {
    int degree, weight;
    degree_weight() : degree(0), weight(0) {}
    degree_weight(int degree, int weight) : degree(degree), weight(weight) {}
    bool operator<(const degree_weight& r) const {
      return degree < r.degree;
    }
    degree_weight& operator+=(const degree_weight& r) {
      degree += r.degree;
      weight += r.weight;
      return *this;
    }
    degree_weight operator+(const degree_weight& r) const {
      degree_weight ret = *this;
      return ret += r;
    }
  };
  void build_depth() {
    depth_.resize(sz(binary_tree_edges_), -1);
    parent_weight_.resize(sz(binary_tree_edges_), make_pair(-2, -2));

    FOR(v, num_vertices_) {
      if (depth_[v] >= 0) continue;

      depth_[v] = 0;
      parent_weight_[v] = make_pair(-1, -1);
      queue<int> q;
      q.push(v);
      while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto& to_weight : binary_tree_edges_[u]) {
          int to, weight; tie(to, weight) = to_weight;
          if (depth_[to] >= 0) continue;
          depth_[to] = depth_[u] + 1;
          parent_weight_[to].first = u;
          parent_weight_[to].second = weight;
          // printf("%d - %d\n",u, to);
          q.push(to);
        }
      }
    }
  }

  void contraction(vector<int>& ancestor, vector<unordered_map<V, int>>& contraction_graph_edges, V u, V v, queue<int>& degree_eq_1) {
    assert(!uf_.is_same(u, v));
    V new_vertex = (V)binary_tree_edges_.size();
    u = uf_.root(u), v = uf_.root(v);
    degree_weight uw = uf_.weight(u), vw = uf_.weight(v);
    uf_.unite(u, v);
    if (u != uf_.root(u)) {
      swap(u, v);
      swap(uw, vw);
    }

    int ua = ancestor[u], va = ancestor[v];
    ancestor[u] = new_vertex;

    // 縮約した結果を元にbinary treeの構築を進める
    binary_tree_edges_.emplace_back(); // new_vertex分の確保
    binary_tree_edges_[new_vertex].emplace_back(ua, uw.weight);
    binary_tree_edges_[ua].emplace_back(new_vertex, uw.weight);
    binary_tree_edges_[new_vertex].emplace_back(va, vw.weight);
    binary_tree_edges_[va].emplace_back(new_vertex, vw.weight);

    //縮約した頂点の辺をまとめる
    auto& uset = contraction_graph_edges[u];
    auto& vset = contraction_graph_edges[v];
    int dec_weight = 0;
    for (const auto& edge : vset) {
      V to; int w; tie(to, w) = edge;
      if (to == u) {
        auto it = uset.find(v);
        assert(it != uset.end());
        assert(it->second == w);
        uset.erase(it);
        dec_weight += w * 2;
      } else {
        uset[to] += w;
        auto& toset = contraction_graph_edges[to];
        auto it = toset.find(v);
        toset.erase(it);
        toset[u] += w;
        if (FLAGS_prune_if_degree_eq_1 && sz(toset) == 1) {
          degree_eq_1.push(to);
        }
      }
    }
    vset.clear();
    uf_.add_weight(u, degree_weight(-2, -dec_weight));
  }

public:
  // g の辺を使い、グラフを作成する(勝手にundirectedとして読み替えている)
  random_contraction2(G& g) :
    num_vertices_(g.num_vertices()),
    uf_(g.num_vertices()),
    binary_tree_edges_(g.num_vertices()) {
    //unordered -> weight = 1なので、shuffleでよい
    //重みがあるならBITで。(stochastic acceptanceはupdateありだと使えない)
    vector<unordered_map<V, int>> contraction_graph_edges(g.num_vertices());
    typename G::edge_list_type initial_edges(g.edge_list());
    std::shuffle(initial_edges.begin(), initial_edges.end(), agl::random);
    vector<int> ancestor(g.num_vertices());

    FOR(i, num_vertices_) ancestor[i] = i;

    for (auto edge : initial_edges) {
      V u = edge.first, v = to(edge.second);
      contraction_graph_edges[u].emplace(v, 1);
      contraction_graph_edges[v].emplace(u, 1);
    }
    FOR(v, num_vertices_) uf_.add_weight(v, degree_weight(sz(contraction_graph_edges[v]), sz(contraction_graph_edges[v])));

    queue<int> degree_eq_1;
    int pruned_degree_eq_1 = 0;
    if (FLAGS_prune_if_degree_eq_1) {
      FOR(v, num_vertices_) {
        if (sz(contraction_graph_edges[v]) == 1) degree_eq_1.push(v);
      }
    }

    for (const auto& edge : initial_edges) {
      if (FLAGS_prune_if_degree_eq_1) {
        // 次数1の頂点があればそれを縮約する
        while (!degree_eq_1.empty()) {
          V v = degree_eq_1.front(); degree_eq_1.pop();
          if (contraction_graph_edges[v].size() == 0) continue;
          V u;
          tie(u, std::ignore) = *contraction_graph_edges[v].begin();
          contraction(ancestor, contraction_graph_edges, u, v, degree_eq_1);
          if (contraction_graph_edges[u].size() == 1) degree_eq_1.push(u);

          pruned_degree_eq_1++;
        }
      }

      V u = edge.first, v = to(edge.second);
      if (uf_.is_same(u, v)) continue;
      contraction(ancestor, contraction_graph_edges, u, v, degree_eq_1);
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
      if (u == -1 || v == -1) return 0; // disconnect
      if (depth_[u] > depth_[v]) swap(u, v);
      ans = min(ans, parent_weight_[v].second);
      v = parent_weight_[v].first;
    }
    return ans;
  }

  void single_source_mincut_dfs(V v, V par, int mn, vector<int>& ans) const {
    if (v < num_vertices_) {
      ans[v] = mn;
    }
    for (auto& to : binary_tree_edges_[v]) {
      if (to.first == par) continue;
      int nmn = min(mn, to.second);
      single_source_mincut_dfs(to.first, v, nmn, ans);
    }
  }


  vector<int> single_source_mincut(V u) const {
    vector<int> ans(num_vertices_, numeric_limits<int>::max());
    single_source_mincut_dfs(u, -1, numeric_limits<int>::max(), ans);
    return ans;
  }

private:
  const int num_vertices_;
  weighted_union_find<degree_weight> uf_;

  vector<vector<pair<V, int>>> binary_tree_edges_;
  vector<pair<V, int>> parent_weight_;
  vector<int> depth_;
};

class min_cut_query2 {
public:
  min_cut_query2(G& g) {
    fprintf(stderr, "initialize mincut-query tree ...\n");
    for (int i = 0; i < FLAGS_solver_iter; i++) {
      if (i % 10 == 0) fprintf(stderr, "trees .. %d/%d\n", i, FLAGS_solver_iter);
      solvers_.emplace_back(g);
    }
    fprintf(stderr, "completed.\n");
  }

  int query(int u, int v) const {
    int ans = numeric_limits<int>::max();
    for (const auto& solver : solvers_) {
      ans = min(ans, solver.query(u, v));
    }
    return ans;
  }

  vector<int> single_source_mincut(V v) const {
    vector<int> ret;
    for (const auto& solver : solvers_) {
      auto cur_ans = solver.single_source_mincut(v);
      if (sz(ret) == 0) ret = move(cur_ans);
      else {
        FOR(i, sz(ret)) ret[i] = min(ret[i], cur_ans[i]);
      }
    }
    return ret;
  }
private:
  vector<random_contraction2> solvers_;
};

void test(G&& g) {
  dinic_twosided dc(g);
  min_cut_query mcq(g);
  // mcq.debug_output_graph();

  JLOG_PUT("argv.mincut_tree_iter", FLAGS_solver_iter);

  xorshift64star gen_node(FLAGS_node_pair_random_seed);
  int unmatch = 0;
  for (int counter = 0; counter < FLAGS_num_query; counter++) {
    V s = gen_node() % g.num_vertices();
    V t = gen_node() % (g.num_vertices() - 1);
    if (s <= t) t++;
    if (counter % 100 == 0) {
      fprintf(stderr, "count/unmatch/all : %d/%d/%d, \n", counter, unmatch, FLAGS_num_query);
    }
    bool is_matched = check_min_cut_query(mcq, dc, s, t, g);
    if (!is_matched) unmatch++;
  }
  JLOG_PUT("result.all", FLAGS_num_query);
  JLOG_PUT("result.match", (FLAGS_num_query - unmatch));
  JLOG_PUT("result.unmatch", unmatch);
}

void test2(G&& g) {
  auto r = agl::random;
  agl::random = r;
  min_cut_query q1(g);
  agl::random = r;
  min_cut_query2 q2(g);

  const int n = g.num_vertices();
  FOR(v, n) for (int s = v + 1; s < n; s++) {
    auto ans1 = q1.query(v, s);
    auto ans2 = q2.query(v, s);
    if (ans1 != ans2) {
      cout << "?";
    }
  }
}

G to_directed_graph(G&& g) {
  vector<pair<V, V>> ret;
  for (auto& e : g.edge_list()) {
    if (e.first < to(e.second)) ret.emplace_back(e.first, to(e.second));
  }
  return G(ret);
}

string graph_name() {
  string x = FLAGS_graph;
  string ret;
  for (int i = sz(x) - 1; i >= 0; i--) {
    if (x[i] == '/' || x[i] == '\\') break;
    ret.push_back(x[i]);
  }
  reverse(ret.begin(), ret.end());
  return ret;
}

DEFINE_string(single_source_mincut_output, "random_contraction_ssm.data", "");
void single_source_mincut(G&& g) {
  min_cut_query mcq(g);
  size_t max_deg = 0;
  int max_deg_v = -1;
  FOR(v, g.num_vertices()) {
    if (g.degree(v, D(0)) + g.degree(D(1)) > max_deg) {
      max_deg = g.degree(v, D(0)) + g.degree(D(1));
      max_deg_v = v;
    }
  }

  FILE* fp = fopen(FLAGS_single_source_mincut_output.c_str(), "w");
  FOR(v, g.num_vertices()) {
    if (v == max_deg_v) continue;
    int xx = mcq.query(max_deg_v, v);
    fprintf(fp, "%d %d %d\n", max_deg_v, v, xx);
  }
  fclose(fp);
}

void single_source_mincut_with_deg(G&& g) {
  size_t max_deg = 0;
  int max_deg_v = -1;
  FOR(v, g.num_vertices()) {
    if (g.degree(v) > max_deg) {
      max_deg = g.degree(v);
      max_deg_v = v;
    }
  }

  FILE* fp = fopen(FLAGS_single_source_mincut_output.c_str(), "w");
  FOR(v, g.num_vertices()) {
    if (v == max_deg_v) continue;
    int xx = g.degree(v, D(0)) + g.degree(v, D(1));
    fprintf(fp, "%d %d %d\n", max_deg_v, v, xx);
  }
  fclose(fp);
}

void tester() {
  test(to_directed_graph(built_in_graph("karate_club")));
  test(to_directed_graph(built_in_graph("dolphin")));
  test(to_directed_graph(built_in_graph("ca_grqc")));

  return;
}

int main(int argc, char **argv) {
  // JLOG_INIT(&argc, argv); called in "easy_cui_init"
  // tester();
  G g = easy_cui_init(argc, argv);
  g = to_directed_graph(std::move(g));

  if (FLAGS_method == "test") {
    test(move(g));
  } else if (FLAGS_method == "single_source_mincut") {
    single_source_mincut(move(g));
  } else if (FLAGS_method == "single_source_mincut_with_deg") {
    single_source_mincut_with_deg(move(g));
  } else {
    fprintf(stderr, "unrecognized option '-method=%s'\n", FLAGS_method.c_str());
    exit(-1);
  }


  return 0;
}
