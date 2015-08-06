#include "base/base.h"
#include "graph/graph.h"
#include "edge_centrality/edge_centrality.h"
#include <gflags/gflags.h>
using namespace std;
using namespace agl;

DEFINE_string(graph, "-", "input graph");

int main(int argc, char **argv) {
  google::ParseCommandLineFlags(&argc, &argv, true);

  G g = read_graph_tsv(FLAGS_graph.c_str());

  auto ebc = edge_betweenness_centrality(g);
  for (auto i : ebc) {
    cout << i.first.first << "\t" << i.first.second << "\t" << i.second << endl;
  }

  return 0;
}
