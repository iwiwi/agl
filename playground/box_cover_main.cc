#include "easy_cui.h"
#include "box_cover.h"
using namespace std;

DEFINE_int32(rad_min, 0, "minimum radius");
DEFINE_int32(rad_max, 10, "maximum radius");
DEFINE_string(method, "sketch", "using method");
DEFINE_double(final_coverage, 0.98, "coverage");
DEFINE_int32(pass, 1, "Number of multi-pass");
DEFINE_int32(sketch_k, 1024, "sketch k");

int main(int argc, char** argv) {
  auto str_replace = [](string& str, const string& from, const string& to) {
    string::size_type pos = 0;
    while (pos = str.find(from, pos), pos != string::npos) {
      str.replace(pos, from.length(), to);
      pos += to.length();
    }
  };

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
  }

  G g(es);
  pretty_print(g);
  JLOG_ADD_OPEN("graph_info") {
    JLOG_PUT("vertices", g.num_vertices());
    JLOG_PUT("edges", g.num_edges());
    JLOG_PUT("graph", FLAGS_graph);
  }

  const char* slash = strrchr(FLAGS_graph.c_str(), '/');
  string graph_name = slash ? slash + 1 : FLAGS_graph;
  str_replace(graph_name, " ", "-");

  string experiment_name = FLAGS_method;
  if (FLAGS_method == "sketch") {
    experiment_name.append("_k=" + to_string(FLAGS_sketch_k));
    experiment_name.append("_pass=" + to_string(FLAGS_pass));
    experiment_name.append("_coverage=" + to_string(FLAGS_final_coverage));

    if (FLAGS_pass == 1) {
      FLAGS_final_coverage = 1.0;
    }
    JLOG_ADD_OPEN("algorithms") {
      JLOG_PUT("name", experiment_name);
      JLOG_PUT("k", to_string(FLAGS_sketch_k));
      JLOG_PUT("pass", to_string(FLAGS_pass));

      for (W rad = FLAGS_rad_min; rad <= FLAGS_rad_max; ++rad) {
        vector<V> res;
        JLOG_ADD_BENCHMARK("time") res = box_cover_sketch(
            g, rad, FLAGS_sketch_k, FLAGS_pass, FLAGS_final_coverage);
        JLOG_ADD("size", res.size());
        JLOG_ADD("coverage", coverage(g, res, rad));
        if (res.size() == 1) break;
      }
    }
  } else if (FLAGS_method == "memb") {
    JLOG_ADD_OPEN("algorithms") {
      JLOG_PUT("name", "MEMB");

      for (W rad = FLAGS_rad_min; rad <= FLAGS_rad_max; ++rad) {
        vector<V> res;
        JLOG_ADD_BENCHMARK("time") res = box_cover_memb(g, rad);
        JLOG_ADD("size", res.size());
        JLOG_ADD("coverage", coverage(g, res, rad));
        if (res.size() == 1) break;
      }
    }
  } else if (FLAGS_method == "burning") {
    JLOG_ADD_OPEN("algorithms") {
      JLOG_PUT("name", "Burning");

      for (W rad = FLAGS_rad_min; rad <= FLAGS_rad_max; ++rad) {
        vector<V> res;
        JLOG_ADD_BENCHMARK("time") res = box_cover_burning(g, rad);
        JLOG_ADD("size", res.size());
        JLOG_ADD("coverage", coverage(g, res, rad));
        if (res.size() == 1) break;
      }
    }
  } else {
    cerr << "Unknown method: " << FLAGS_method << endl;
    return 0;
  }

  JLOG_INSERT_FILENAME("_rad_max=" + to_string(FLAGS_rad_max));
  JLOG_INSERT_FILENAME("_" + graph_name);
  JLOG_INSERT_FILENAME(experiment_name);
}