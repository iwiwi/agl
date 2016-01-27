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
  ofstream os(FLAGS_validation_data_path.c_str(), ios_base::out);
  gf.print_gomory_hu_tree(os);
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
  Gusfield3 gf3(g);
  Gusfield4 gf4(g);
  auto l = f(gf3);
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
      T gf(g);
      JLOG_ADD("getcap_counter", getcap_counter);
      JLOG_ADD("addcap_counter", addcap_counter);
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
