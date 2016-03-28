#include "easy_cui.h"
#include "shortest_path/dynamic_pruned_landmark_labeling.h"

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);
  dynamic_pruned_landmark_labeling<16> dpll;
  V num_v = g.num_vertices();
  std::vector<bool> used(num_v, false);

  JLOG_PUT_BENCHMARK("result.load") { dpll.load_graph(g); }
  JLOG_PUT_BENCHMARK("result.bit_parallel") {
    dpll.bit_parallel_bfs(used, g.num_edges());
  }
  JLOG_PUT_BENCHMARK("result.pruned_bfs") {
    for (V root = 0; root < num_v; ++root) {
      if (used[root]) continue;
      dpll.pruned_bfs(root, 0, used);
      dpll.pruned_bfs(root, 1, used);
      used[root] = true;
    }
  }

  JLOG_PUT("result.average_label_size", (dpll.total_label_num() / num_v));
  JLOG_PUT_BENCHMARK("result.query_time_1M") {
    for (size_t i = 0; i < 1000 * 1000; ++i) {
      dpll.query_distance(g, agl::random(num_v), agl::random(num_v));
    }
  }
  JLOG_PUT_BENCHMARK("result.add_time_1K") {
    for (size_t i = 0; i < 1000; ++i) {
      dpll.add_edge(g, agl::random(num_v), agl::random(num_v));
    }
  }

  return 0;
}
