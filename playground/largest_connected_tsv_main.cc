#include "easy_cui.h"
#include "box_cover.h"
using namespace std;

DEFINE_int32(rad_min, 1, "minimum radius");
DEFINE_int32(rad_max, 10, "maximum radius");
DEFINE_string(method, "sketch", "using method");
DEFINE_double(final_coverage, 1.0, "coverage");
DEFINE_int32(pass, 1, "Number of multi-pass");
DEFINE_int32(sketch_k, 128, "sketch k");
DEFINE_bool(rad_analytical, false, "Using analytical diameters for rads");
DEFINE_double(upper_param, 1.0, "size_upper_bound=upper_param*n*k");
DEFINE_string(exp_tag, "", "experiment name");

unweighted_edge_list extract_maximal_connected(const G& g_pre) {
  unweighted_edge_list es;
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

  vector<V>& max_c = connected_sets[max_g];
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
  return es;
}

int main(int argc, char** argv) {
  unweighted_edge_list es;
  {
    G g_pre = easy_cui_init(argc, argv);
    CHECK_MSG(FLAGS_force_undirected, "undirected only!!!");
    es = extract_maximal_connected(g_pre);
  }
  es = make_undirected(es);
  for (const auto& p : es) cout << p.first << "\t" << p.second << endl;
}
