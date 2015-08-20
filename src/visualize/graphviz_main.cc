#include "easy_cui.h"
#include "graphviz.h"
DEFINE_string(out, "out", "output prefix");
DEFINE_string(command, "dot", "graphviz command");

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);
  graphviz(g, FLAGS_out.c_str(), FLAGS_command.c_str());
  return 0;
}
