#include <easy_cui.h>

#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())

#include "dinic_twosided.h"
#include "dinic_onesided.h"
#include "ConnectedComponentsFilter.h"

DEFINE_string(method, "oneside", "oneside, twoside");
DEFINE_int32(num_query, 1000, "");
DEFINE_int64(node_pair_random_seed, 922337203685477583LL, "");

// 次数がそこそこ低い頂点同士で、dinicの速度を測る
template<class maxflow_T>
void test(G& g){
  xorshift64star gen_node(FLAGS_node_pair_random_seed);
  
  maxflow_T dc(g);
  vector<V> vs;
  for(V v : make_irange(g.num_vertices())) {
    vs.push_back(v);
  }
  for (int counter = 0; counter < FLAGS_num_query; counter++) {
      int sid = gen_node() % sz(vs);
      int tid = gen_node() % (sz(vs) - 1);
      if (sid <= tid) tid++;
      V s = vs[sid] , t = vs[tid];
      if (counter % 100 == 0) {
        fprintf(stderr, "count/all : %d/%d, \n", counter, FLAGS_num_query);
      }
      dc.reset_graph();
      int flow = dc.max_flow(s, t);
      fprintf(stderr, "(%d,%d) = %d\n", s, t, flow);
    }

}

G to_directed_graph(G g) {
  vector<pair<V, V>> ret;
  for (auto& e : g.edge_list()) {
    if (e.first < to(e.second)) ret.emplace_back(e.first, to(e.second));
  }
  return G(ret);
}

int main(int argc, char** argv) {

  G g = easy_cui_init(argc, argv);
  // g = to_directed_graph(g);

  if (FLAGS_method == "oneside") {
    test<dinic_onesided>(g);
  } else if (FLAGS_method == "twoside") {
    test<dinic_twosided>(g);
  } else {
    fprintf(stderr, "unrecognized option '%s'\n", FLAGS_method.c_str());
    exit(-1);
  }
}
