#include "easy_cui.h"

#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())

#include "dinic_twosided.h"
#include "dinic_naive.h"
#include "random_contraction.h"

DEFINE_int32(num_query, 1000, "");
DEFINE_int64(node_pair_random_seed, 922337203685477583LL, "");
DEFINE_string(method, "test", "test, single_source_mincut");

bool check_min_cut_query(min_cut_query2& mcq2, dinic_twosided& dc, int s, int t, G& g) {
  dc.reset_graph();
  int naive_w = dc.max_flow(s, t);
  int mcq_w = mcq2.query(s, t);
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

void test(G&& g) {
  min_cut_query2 mcq2(g);
  dinic_twosided dc(g);
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
    bool is_matched = check_min_cut_query(mcq2, dc, s, t, g);
    if (!is_matched) unmatch++;
  }
  JLOG_PUT("result.all", FLAGS_num_query);
  JLOG_PUT("result.match", (FLAGS_num_query - unmatch));
  JLOG_PUT("result.unmatch", unmatch);
}

G to_directed_graph(G&& g) {
  vector<pair<V, V>> ret;
  for (auto& e : g.edge_list()) {
    if (e.first < to(e.second)) ret.emplace_back(e.first, to(e.second));
    else if (to(e.second) < e.first) ret.emplace_back(to(e.second), e.first);
  }
  sort(ret.begin(), ret.end());
  ret.erase(unique(ret.begin(), ret.end()), ret.end());
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
  min_cut_query2 mcq(g);
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
    if (v % 10000 == 0) cerr << v << endl;
    auto ans = mcq.query(max_deg_v, v);
    fprintf(fp, "%d %d %d\n", max_deg_v, v, ans);
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
  exit(0);
}

int main(int argc, char **argv) {
  // tester();
  // JLOG_INIT(&argc, argv); called in "easy_cui_init"
  G g = easy_cui_init(argc, argv);
  g = to_directed_graph(move(g));
  pretty_print(g);

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
