#include "easy_cui.h"
#include "graphviz.h"
#include "edge_centrality/edge_centrality.h"
DEFINE_string(out, "out.png", "output PNG image");

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);
  edge_centrality_map ec = edge_betweenness_centrality_naive(g);
  graphviz_draw_edge_centrality(g, ec, FLAGS_out.c_str());
  return 0;
}
