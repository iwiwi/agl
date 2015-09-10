#include "vertex_centrality.h"
using namespace std;

namespace agl {
vector<double> vertex_centrality_degree(const G &g) {
  vector<double> vc(g.num_vertices());
  for (V v : g.vertices()) {
    vc[v] = g.degree(v);
  }
  return vc;
}
}  // namespace agl
