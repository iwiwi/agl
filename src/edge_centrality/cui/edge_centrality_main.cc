#include "easy_cui.h"
#include "edge_centrality/edge_centrality.h"

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);

  auto ebc = edge_betweenness_centrality(g);
  if (FLAGS_force_undirected) {
    ebc = merge_edge_centrality_map_entries_for_undirected_graph(ebc);
  }

  for (auto i : ebc) {
    cout << i.first.first << "\t" << i.first.second << "\t" << i.second << endl;
  }

  return 0;
}
