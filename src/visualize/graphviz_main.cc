#include "base/base.h"
#include "graph/graph.h"
#include "graphviz.h"
#include <gflags/gflags.h>
using namespace std;
using namespace agl;

DEFINE_string(graph, "-", "input graph");
DEFINE_string(out, "out", "output prefix");
DEFINE_string(command, "dot", "graphviz command");

int main(int argc, char **argv) {
  google::ParseCommandLineFlags(&argc, &argv, true);

  G g = read_graph_tsv(FLAGS_graph.c_str());
  graphviz(g, FLAGS_out.c_str(), FLAGS_command.c_str());

  return 0;
}
