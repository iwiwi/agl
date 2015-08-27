#include "easy_cui.h"
#include "graphviz.h"
DEFINE_string(out, "out.png", "output prefix");

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);
  graphviz_draw_graph(g, FLAGS_out.c_str());
  return 0;
}
