#include <easy_cui.h>

#ifdef _WIN32
#pragma comment(linker, "/STACK:3200000") 
#endif

DEFINE_string(method, "gusfield", "test, gusfield, print_gomory_hu_tree");

#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())

#include "bi_dinitz.h"
#include "ConnectedComponentsFilter.h"
#include "TwoEdgeCCFilter.h"
#include "OptimizedGusfieldWith2ECC.h"
#include "OptimizedGusfieldWith2ECC_slow.h"
#include "PlainGusfield.h"
#include "PlainGusfield_bi_dinitz.h"

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

using Gusfield = TwoEdgeCCFilter<OptimizedGusfieldWith2ECC>;
using Gusfield_slow = TwoEdgeCCFilter<OptimizedGusfieldWith2ECC_slow>;
DEFINE_string(gomory_hu_builder, "Gusfield", "Gusfield, Gusfield_slow");

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

DEFINE_string(validation_data_path, "", "");

template<class gomory_hu_tree_t>
void print_gomory_hu_tree(G&& g) {

  fprintf(stderr, "print_gomory_hu_tree : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);

  if (FLAGS_validation_data_path == "") {
    auto gname = graph_name();
    FLAGS_validation_data_path = gname + ".tree";
  }
  gomory_hu_tree_t* gf;
  JLOG_PUT_BENCHMARK("gusfield_time") {
    gf = new gomory_hu_tree_t(g);
  }

  ofstream os(FLAGS_validation_data_path.c_str(), ios_base::out);
  gf->print_gomory_hu_tree(os);
  delete gf;
}

map<int, vector<pair<V, V>>> load(const string& path) {
  ifstream is(path.c_str());
  map<int, vector<pair<V, V>>> ret;
  int s, v, weight;
  while (is >> s >> v >> weight) {
    ret[weight].emplace_back(s, v);
  }
  return ret;
}

template<class T>
map<int, vector<pair<V, V>>> f(T& gf) {
  stringstream ss;
  gf.print_gomory_hu_tree(ss);
  map<int, vector<pair<V, V>>> ret;
  int s, v, weight;
  while (ss >> s >> v >> weight) {
    ret[weight].emplace_back(s, v);
  }
  return ret;
}

void test(G&& g) {
  Gusfield gf(g);
  Gusfield_slow gf4(g);
  auto l = f(gf);
  auto r = f(gf4);

  using ull = unsigned long long;
  const int n = g.num_vertices();
  vector<ull> zobrist_hash(n);
  FOR(i, n) zobrist_hash[i] = (ull(agl::random()) << 32) | ull(agl::random());

  union_find ufl(n), ufr(n);
  vector<ull> hl = zobrist_hash, hr = zobrist_hash;

  vector<int> weight;
  for (auto& kv : l) weight.push_back(kv.first);
  reverse(weight.begin(), weight.end());

  auto on_mismatched = [&](V u) {
    set<int> onlyl, onlyr, intersect;
    FOR(a, n) if (ufl.is_same(u, a)) onlyl.insert(a);
    FOR(a, n) if (ufr.is_same(u, a)) {
      if (onlyl.count(a)) onlyl.erase(a), intersect.insert(a);
      else onlyr.insert(a);
    }
    int cnt;
    printf("left: ");
    cnt = 0;
    for (auto x : onlyl) {
      printf("%d, ", x);
      cnt++;
      if (cnt >= 10) break;
    }
    puts("");
    printf("right: ");
    cnt = 0;
    for (auto x : onlyr) {
      printf("%d, ", x);
      cnt++;
      if (cnt >= 10) break;
    }
    puts("");
    printf("intersect: ");
    cnt = 0;
    for (auto x : intersect) {
      printf("%d, ", x);
      cnt++;
      if (cnt >= 10) break;
    }
    puts("");
  };

  for (const int w : weight) {
    if (sz(l[w]) != sz(r[w])) printf("*");
    printf("w = %d, L : %d , R : %d", w, sz(l[w]), sz(r[w]));
    puts("");
  }

  for (const int w : weight) {
    // CHECK(sz(l[w]) == sz(r[w]));
    for (auto& uv : l[w]) {
      int u, v; tie(u, v) = uv;
      u = ufl.root(u), v = ufl.root(v);
      CHECK(u != v);
      ufl.unite(u, v);
      ull new_hash = hl[u] ^ hl[v];
      hl[ufl.root(u)] = new_hash;
    }

    for (auto& uv : r[w]) {
      int u, v; tie(u, v) = uv;
      u = ufr.root(u), v = ufr.root(v);
      CHECK(u != v);
      ufr.unite(u, v);
      ull new_hash = hr[u] ^ hr[v];
      hr[ufr.root(u)] = new_hash;
    }

    for (auto& uv : l[w]) {
      int u, v; tie(u, v) = uv;
      if (hl[ufl.root(u)] != hr[ufr.root(u)]) {
        on_mismatched(u);
        exit(-1);
      }
      if (hl[ufl.root(v)] != hr[ufr.root(v)]) {
        on_mismatched(v);
        exit(-1);
      }
    }
    for (auto& uv : r[w]) {
      int u, v; tie(u, v) = uv;
      if (hl[ufl.root(u)] != hr[ufr.root(u)]) {
        on_mismatched(u);
        exit(-1);
      }
      if (hl[ufl.root(v)] != hr[ufr.root(v)]) {
        on_mismatched(v);
        exit(-1);
      }
    }
  }
}

void tester() {
  test(to_directed_graph(built_in_graph("karate_club")));
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
      T gf(g);
    }
  } else if (FLAGS_method == "print_gomory_hu_tree") {
    print_gomory_hu_tree<T>(std::move(g));
  } else {
    fprintf(stderr, "unrecognized option '%s'\n", FLAGS_method.c_str());
    exit(-1);
  }

  JLOG_PUT("try_greedy_tree_packing", FLAGS_try_greedy_tree_packing);
  JLOG_PUT("getcap_counter", getcap_counter);
  JLOG_PUT("addcap_counter", addcap_counter);
  JLOG_PUT("preflow_eq_degree", preflow_eq_degree);
  JLOG_PUT("flow_eq_0", flow_eq_0);
  JLOG_PUT("gtp_edge_count_all", gtp_edge_count);
  JLOG_PUT("gtp_edge_miss_all", gtp_edge_miss);
  JLOG_PUT("gtp_edge_use_all", gtp_edge_use);
}

DEFINE_bool(to_directed_graph,true, "");
DEFINE_string(write_directed_graph_name,"", "");

int main(int argc, char** argv) {

  // tester();

  G g = easy_cui_init(argc, argv);
  fprintf(stderr, "easy_cui_init : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);
  if(FLAGS_graph.find(".directed") == string::npos) {
    g = to_directed_graph(std::move(g));
    fprintf(stderr, "load graph : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);
  }

  if(FLAGS_gomory_hu_builder == "write_directed_graph") {
    string output = FLAGS_write_directed_graph_name;
    if(output == "") {
      output = FLAGS_graph + ".directed"; 
    }
    write_graph_binary(g, output.c_str());
    exit(0);
  }

  if (FLAGS_gomory_hu_builder == "PlainGusfield") { 
    main_<PlainGusfield>(std::move(g));
  } else if (FLAGS_gomory_hu_builder == "PlainGusfield_bi_dinitz") { 
    main_<PlainGusfield_bi_dinitz>(std::move(g));
  } else if (FLAGS_gomory_hu_builder == "Gusfield") {
    main_<Gusfield>(std::move(g));
  } else if (FLAGS_gomory_hu_builder == "Gusfield_slow") {
    main_<Gusfield_slow>(std::move(g));
  } else {
    fprintf(stderr, "unrecognized option -gomory_hu_builder='%s'\n", FLAGS_gomory_hu_builder.c_str());
    exit(-1);
  }
}
