#include "dijkstra.h"
using namespace std;

namespace agl {
vector<pair<W, double>> single_source_distance_with_num_paths(const G &g, V s, D dir) {
  vector<pair<W, double>> dps(g.num_vertices(), {kInfW, 0});
  auto h = make_dijkstra_heap(g);

  h.decrease(s, 0);
  dps[s] = make_pair(0, 1);
  while (!h.empty()) {
    V v = h.top_vertex();
    W d = h.top_weight();
    h.pop();

    for (auto e : g.edges(v, dir)) {
      V tv = to(e);
      W td = d + weight(e);
      h.decrease(tv, td);
      if (is_lt(td, dps[tv].first)) dps[tv] = make_pair(td, dps[v].second);
      else if (is_eq(td, dps[tv].first)) dps[tv].second += dps[v].second;
    }
  }

  return dps;
}
}  // namespace agl
