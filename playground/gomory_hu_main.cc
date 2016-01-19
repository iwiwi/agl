#include <easy_cui.h>

#ifdef _WIN32
#pragma comment(linker, "/STACK:3200000") 
#endif

DEFINE_int32(num_query, 1000, "");
DEFINE_int64(node_pair_random_seed, 922337203685477583LL, "");
DEFINE_string(method, "gusfield", "test, gusfield, print_gomory_hu_tree");

#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())

#include "dinic_naive.h"
#include "dinic_twosided.h"
#include "ConnectedComponentsFilter.h"
#include "TwoEdgeCCFilter.h"
#include "OptimizedGusfield.h"
#include "OptimizedGusfieldWith2ECC.h"
#include "OptimizedGusfieldWith2ECC2.h"

G to_directed_graph(G g) {
  vector<pair<V, V>> ret;
  for (auto& e : g.edge_list()) {
    if (e.first < to(e.second)) ret.emplace_back(e.first, to(e.second));
  }
  return G(ret);
}

using Gusfield3 = TwoEdgeCCFilter<OptimizedGusfieldWith2ECC>;
using Gusfield4 = TwoEdgeCCFilter<OptimizedGusfieldWith2ECC2>;
DEFINE_string(gomory_fu_builder, "Gusfield4", "Gusfield3, Gusfield4");

void aggregate_weight(const G& g) {
  Gusfield3 gf3(g);
  gf3.aggregate_gomory_hu_tree_weight();
}

DEFINE_string(validation_data_path, "validate.data", "");

template<class gomory_hu_tree_t>
void print_gomory_hu_tree(G&& g) {
  gomory_hu_tree_t gf(g);
  FILE* fp = fopen(FLAGS_validation_data_path.c_str(), "w");
  gf.print_gomory_hu_tree(fp);
  fclose(fp);
}

void test(G&& g) {
  Gusfield3 gf3(g);
  Gusfield4 gf4(g);

  int counter = 0;
  const int all = g.num_vertices() * (g.num_vertices() - 1) / 2;
  FOR(s, g.num_vertices()) for (int t = s + 1; t < g.num_vertices(); t++) {
    if (counter % 100 == 0) {
      fprintf(stderr, "count/all = %d/%d\n", counter, all);
    }
    int flow3 = gf3.query(s, t);
    int flow4 = gf4.query(s, t);
    printf("(%d,%d) : gf3 = %d, gf4 = %d\n", s, t, flow3, flow4);
    if (flow3 != flow4) {
      puts("?");
      int x; cin >> x;
    }
    counter++;
  }
}

void tester() {
  FLAGS_validation_data_path = "Gusfield3.data";
  print_gomory_hu_tree<Gusfield3>(to_directed_graph(built_in_graph("karate_club")));
  FLAGS_validation_data_path = "Gusfield4.data";
  print_gomory_hu_tree<Gusfield4>(to_directed_graph(built_in_graph("karate_club")));
  test(to_directed_graph(built_in_graph("dolphin")));
  test(to_directed_graph(built_in_graph("ca_grqc")));
  exit(0);
}

template<class T>
void main_(G&& g) {
  if (FLAGS_method == "test") {
    test(std::move(g));
  } else if (FLAGS_method == "gusfield") {
    JLOG_PUT_BENCHMARK("gusfield_time") {
      Gusfield3 gf(g);
    }
  } else if (FLAGS_method == "print_gomory_hu_tree") {
    print_gomory_hu_tree<T>(std::move(g));
  } else {
    fprintf(stderr, "unrecognized option '%s'\n", FLAGS_method.c_str());
    exit(-1);
  }
}

int main(int argc, char** argv) {

  // tester();

  G g = easy_cui_init(argc, argv);
  g = to_directed_graph(g);

  if (FLAGS_gomory_fu_builder == "Gusfield3") {
    main_<Gusfield3>(std::move(g));
  } else if (FLAGS_gomory_fu_builder == "Gusfield4") {
    main_<Gusfield4>(std::move(g));
  } else {
    fprintf(stderr, "unrecognized option -gomory_fu_builder='%s'\n", FLAGS_gomory_fu_builder.c_str());
    exit(-1);
  }
}
