#include "vertex_centrality.h"
#include "shortest_path/shortest_path.h"
using namespace std;

namespace agl {
vector<double> vertex_centrality_degree(const G &g) {
  vector<double> vc(g.num_vertices());
  for (V v : g.vertices()) {
    vc[v] = g.degree(v);
  }
  return vc;
}

vector<double> vertex_centrality_classic_closeness(const G& g) {
  vector<double> cc(g.num_vertices());
  for (V v : g.vertices()) {
    auto ds = single_source_distance(g, v);
    for (V u : g.vertices()) cc[u] += ds[u];
  }
  for (auto &x : cc) x = 1 / x;
  return cc;
}
}  // namespace agl
