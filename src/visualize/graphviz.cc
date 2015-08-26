#include "graphviz.h"
using namespace std;

DEFINE_string(graphviz_engine, "dot", "dot (default), fdp, neato, twopi, circo, sfdp, ...");

namespace agl {
void graphviz(const G &g, const char *filename, const char *command) {
  string dot_filename = string(filename) + ".dot";
  {
    ofstream ofs(dot_filename.c_str());
    ofs << "digraph G {";
    for (V v : g.vertices()) {
      for (E e : g.edges(v)) {
        ofs << "  " << v << " -> " << to(e) << ";" << endl;
      }
    }
    ofs << "}";
  }
  {
    string cmd = command + string(" ") + dot_filename + " -T png -o " + filename + ".png";
    system(cmd.c_str());
  }
}
}  // namespace agl
