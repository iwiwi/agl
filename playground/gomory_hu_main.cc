#include <easy_cui.h>

#ifdef _WIN32
#pragma comment(linker, "/STACK:3200000") 
#endif

#ifdef _WIN32
template<typename GraphType = G>
GraphType easy_cui_init(int argc, char **argv) {
  JLOG_INIT(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_type == "auto")
    FLAGS_type = guess_type();

  //agl形式は隣接リストを介すると効率が悪いため、直接GraphTypeを生成する
  if (FLAGS_type == "agl") {
    auto g = read_graph_binary<G>(FLAGS_graph.c_str());
    if (FLAGS_force_undirected) {
      g = GraphType(make_undirected(g.edge_list()));
    }
    pretty_print(g);
    return g;
  }

  G::edge_list_type es;
  if (FLAGS_type == "tsv") {
    es = read_edge_list_tsv(FLAGS_graph.c_str());
  } else if (FLAGS_type == "built_in") {
    es = built_in_edge_list(FLAGS_graph.c_str());
  } else if (FLAGS_type == "gen") {
    istringstream iss(FLAGS_graph);
    string family;
    iss >> family;
    if (family == "barbell") {
      V n;
      if (!(iss >> n)) n = 4;
      es = generate_barbell(n);
    } else if (family == "grid") {
      size_t r, c;
      if (!(iss >> r)) r = 4;
      if (!(iss >> c)) c = r;
      es = generate_grid(r, c);
    } else if (family == "erdos_renyi") {
      V n;
      double d;
      if (!(iss >> n)) n = 10;
      if (!(iss >> d)) d = 3.0;
      es = generate_erdos_renyi(n, d);
    } else if (family == "random_planar") {
      V n;
      size_t e;
      if (!(iss >> n)) n = 10;
      if (!(iss >> e)) e = 25;
      es = generate_random_planar(n, e);
    } else if (family == "cycle") {
      V n;
      if (!(iss >> n)) n = 10;
      es = generate_cycle(n);
    } else if (family == "ba") {
      V n, m;
      if (!(iss >> n)) n = 10;
      if (!(iss >> m)) m = 5;
      es = generate_ba(n, m);
    } else if (family == "dms") {
      V n, m, k0;
      if (!(iss >> n)) n = 10;
      if (!(iss >> m)) m = 5;
      if (!(iss >> k0)) k0 = -2;
      es = generate_dms(n, m, k0);
    } else if (family == "hk") {
      V n, m;
      double p;
      if (!(iss >> n)) n = 10;
      if (!(iss >> m)) m = 5;
      if (!(iss >> p)) p = 0.5;
      es = generate_hk(n, m, p);
    } else if (family == "ws") {
      V n, d;
      double p;
      if (!(iss >> n)) n = 10;
      if (!(iss >> d)) d = 4;
      if (!(iss >> p)) p = 0.5;
      es = generate_ws(n, d, p);
    } else if (family == "kronecker") {
      int scale, n;
      size_t avg_deg;
      if (!(iss >> scale)) scale = 5;
      if (!(iss >> avg_deg)) avg_deg = 16;
      vector<vector<double>> mat;
      if (!(iss >> n)) {
        n = 2;
        mat = vector<vector<double>>(n, vector<double>(n));
        mat[0][0] = 0.57;
        mat[0][1] = 0.19;
        mat[1][0] = 0.19;
        mat[1][1] = 0.05;
      }
      for (int i = 0; i < n * n; ++i) {
        if (!(iss >> mat[i / n][i % n])) mat[i / n][i % n] = 1.0 / (n * n);
      }
      es = generate_kronecker(scale, avg_deg, mat);
    } else if (family == "flower") {
      V required, u, v;
      if (!(iss >> required)) required = 44;
      if (!(iss >> u)) u = 2;
      if (!(iss >> v)) v = 2;
      es = generate_uv_flower(required, u, v);
    } else if (family == "shm") {
      V required_num, initial_num;
      int t;
      if (!(iss >> required_num)) required_num = 101;
      if (!(iss >> initial_num)) initial_num = 5;
      if (!(iss >> t)) t = 2;
      es = generate_shm(required_num, initial_num, t);
    } else {
      FAIL_MSG("Unknown generator family: " + family);
    }
  }

  if (FLAGS_force_undirected) es = make_undirected(es);

  GraphType g(es);
  pretty_print(g);
  return g;
}
#endif

DEFINE_int32(num_query, 1000, "");
DEFINE_int64(node_pair_random_seed, 922337203685477583LL, "");
DEFINE_string(method, "gusfield", "both, gomory_hu, gusfield");

#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())

#include "dinic_naive.h"
#include "dinic_twosided.h"
#include "ConnectedComponentsFilter.h"
#include "OptimizedGusfield.h"
#include "OptimizedGusfieldWith2ECC.h"

G to_directed_graph(G g) {
  vector<pair<V, V>> ret;
  for (auto& e : g.edge_list()) {
    if (e.first < to(e.second)) ret.emplace_back(e.first, to(e.second));
  }
  return G(ret);
}

using Gusfield2 = ConnectedComponentsFilter<OptimizedGusfield<dinic_twosided>>;
using Gusfield3 = OptimizedGusfieldWith2ECC;

bool check_min_cut_query(Gusfield2& gf2, Gusfield3& gf3, int s, int t, const G& g) {
  // int naive_w = dc.max_flow(s, t);
  int gf2_w = gf2.query(s, t);
  int gf3_w = gf3.query(s, t);

  if (gf2_w != gf3_w) {
    fprintf(stderr, "unmatched. (S,T) = (%d,%d), naive = None, gf2_w = %d, gf3_w = %d\n", s, t, gf2_w, gf3_w);
  }
  return gf2_w == gf3_w;
}

void test(const G& g) {
  Gusfield3 gf3(g);
  Gusfield2 gf2(g);
  xorshift64star gen_node(FLAGS_node_pair_random_seed);

  if (g.num_vertices() > 5000) {
    int unmatch = 0;
    for (int counter = 0; counter < FLAGS_num_query; counter++) {
      V s = gen_node() % g.num_vertices();
      V t = gen_node() % (g.num_vertices() - 1);
      if (s <= t) t++;
      if (counter % 100 == 0) {
        fprintf(stderr, "count/unmatch/all : %d/%d/%d, \n", counter, unmatch, FLAGS_num_query);
      }
      bool is_matched = check_min_cut_query(gf2, gf3, s, t, g);
      if (!is_matched) unmatch++;
    }
    JLOG_PUT("result.all", FLAGS_num_query);
    JLOG_PUT("result.match", (FLAGS_num_query - unmatch));
    JLOG_PUT("result.unmatch", unmatch);
  } else {
    int counter = 0;
    int unmatch = 0;
    FOR(s, g.num_vertices()) for (int t = s + 1; t < g.num_vertices(); t++) {
      if (counter % 100 == 0) {
        fprintf(stderr, "count/unmatch/all : %d/%d/%d, \n", counter, unmatch, -1);
      }
      bool is_matched = check_min_cut_query(gf2, gf3, s, t, g);
      if (!is_matched) unmatch++;
      counter++;
    }

    JLOG_PUT("result.all", counter);
    JLOG_PUT("result.match", (counter - unmatch));
    JLOG_PUT("result.unmatch", unmatch);
  }
}

DEFINE_string(validation_data_path, "validate.data", "");

void flow_all_ST_pair(G& g) {
  Gusfield3 gf(g);
  
  FILE* fp = fopen(FLAGS_validation_data_path.c_str(), "w");
  int counter = 0;
  const int all = g.num_vertices() * (g.num_vertices() - 1) / 2;
  FOR(s, g.num_vertices()) for (int t = s + 1; t < g.num_vertices(); t++) {
    if (counter % 100 == 0) {
      fprintf(stderr, "count/all = %d/%d\n", counter, all);
    }
    int flow = gf.query(s, t);
    fprintf(fp, "%d %d %d\n", s, t, flow);
    counter++;
  }
  fclose(fp);
}

void tester() {
  test(to_directed_graph(built_in_graph("ca_grqc")));
  exit(0);
}


int main(int argc, char** argv) {

  G g = easy_cui_init(argc, argv);
  g = to_directed_graph(g);

  if (FLAGS_method == "test") {
    test(g);
  } else if (FLAGS_method == "gusfield") {
    Gusfield3 gf(g);
  } else if(FLAGS_method == "all_pair") {
    flow_all_ST_pair(g);
  } else {
    fprintf(stderr, "unrecognized option '%s'\n", FLAGS_method.c_str());
    exit(-1);
  }
}
