#include "easy_cui.h"
#include "edge_centrality/edge_centrality.h"

DEFINE_int32(n, 10, "n");
DEFINE_double(m, 30, "m");

int main(int argc, char **argv) {
  JLOG_INIT(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);

  G g;
  for (;;) {
    auto el = force_undirected(gen_erdos_renyi(FLAGS_n, FLAGS_m / FLAGS_n));
    g.assign(el);
    if (is_connected(g)) break;
  }

  auto ebc = edge_betweenness_centrality(g);
  for (auto i : ebc) {
    cout << i.first.first << "\t" << i.first.second << "\t" << i.second << endl;
  }

  return 0;
}
