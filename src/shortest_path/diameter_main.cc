#include "easy_cui.h"
#include "diameter.h"

int main(int argc, char** argv) {
  G g = easy_cui_init(argc, argv);
  V d;

  JLOG_PUT_BENCHMARK("diameter") {
    d = diameter(g);
  }

  JLOG_OPEN("graph_info") {
    JLOG_PUT("vertices", g.num_vertices());
    JLOG_PUT("edges", g.num_edges());
    JLOG_PUT("diameter", d);
  }

  return 0;
}
