#include <easy_cui.h>

#ifdef _WIN32
#pragma comment(linker, "/STACK:3200000") 
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

void test(const G& g) {
  Gusfield3 gf3(g);
  gf3.aggregate_gomory_hu_tree_weight();
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
    JLOG_PUT_BENCHMARK("gusfield_time") {
      Gusfield3 gf(g);
    }
  } else if(FLAGS_method == "all_pair") {
    flow_all_ST_pair(g);
  } else {
    fprintf(stderr, "unrecognized option '%s'\n", FLAGS_method.c_str());
    exit(-1);
  }
}
