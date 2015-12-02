#include "easy_cui.h"
#include "box_cover.h"
using namespace std;

DEFINE_int32(rad_min, 0, "minimum radius");
DEFINE_int32(rad_max, 10, "maximum radius");

namespace {
// What portion of the vertices are really covered by the set?
double coverage(const G &g, vector<V> s, W rad) {
  vector<V> covered;
  vector<bool> vis(g.num_vertices());
  queue<pair<V, W>> que;
  for (V center : s) {
    que.push(make_pair(center, 0));
    covered.push_back(center);
    vis[center] = true;
  }
  while (!que.empty()) {
    V v = que.front().first;
    W dist = que.front().second;
    que.pop();
    if (dist == rad) break;
    for (V u : g.neighbors(v)) {
      if (vis[u]) continue;
      vis[u] = true;
      que.push(make_pair(u, dist + 1));
      covered.push_back(u);
    }
  }
  assert(covered.size() <= (size_t)g.num_vertices());

  double coverage = covered.size();
  return coverage / g.num_vertices();
}
}

int main(int argc, char **argv) {
  //
  // Extract maximal connected subgraph & Compress coordinates
  //
  unweighted_edge_list es;
  {
    G g_pre = easy_cui_init(argc, argv);
    CHECK_MSG(FLAGS_force_undirected, "undirected only!!!");

    V num_v = g_pre.num_vertices();
    size_t max_num_v = 0;
    int max_g = 0;
    vector<bool> vis(num_v, false);
    vector<vector<V>> connected_sets;
    for (V v = 0; v < num_v; ++v) {
      if (vis[v]) continue;
      vector<V> connected_v;
      queue<V> que;
      que.push(v);
      vis[v] = true;
      connected_v.push_back(v);
      while (!que.empty()) {
        V u = que.front();
        que.pop();
        for (V x : g_pre.neighbors(u)) {
          if (!vis[x]) {
            vis[x] = true;
            que.push(x);
            connected_v.push_back(x);
          }
        }
      }
      if (max_num_v < connected_v.size()) {
        max_num_v = connected_v.size();
        max_g = connected_sets.size();
        sort(connected_v.begin(), connected_v.end());
        connected_sets.push_back(connected_v);
      }
    }

    vector<V> &max_c = connected_sets[max_g];
    vector<V> inv(num_v, -1);
    for (V v : max_c) {
      inv[v] = (lower_bound(max_c.begin(), max_c.end(), v) - max_c.begin());
    }

    for (pair<V, V> e : g_pre.edge_list()) {
      V from = e.first;
      V to = e.second;
      if (inv[from] == -1) continue;
      es.emplace_back(inv[from], inv[to]);
    }
  }

  G g(es);
  pretty_print(g);

  JLOG_ADD_OPEN("graph_info") {
    JLOG_PUT("vertices", g.num_vertices());
    JLOG_PUT("edges", g.num_edges());
    JLOG_PUT("type", FLAGS_graph);
    JLOG_PUT("graph", FLAGS_graph);
  }

  vector<pair<string, function<vector<V>(const G &, W)>>> algos{
      // {"MEMB", box_cover_memb},
      {"Schneider", box_cover_burning},
      {"Sketch", box_cover_sketch},
  };

  for (auto a : algos) {
    JLOG_ADD_OPEN("algorithms") {
      JLOG_PUT("name", a.first);
      auto f = a.second;

      for (W rad = FLAGS_rad_min; rad <= FLAGS_rad_max; ++rad) {
        vector<V> res;
        JLOG_ADD_BENCHMARK("time") res = f(g, rad);
        JLOG_ADD("size", res.size());
        JLOG_PUT("coverage", coverage(g, res, rad));
      }
    }
  }
}
