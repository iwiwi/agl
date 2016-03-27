#include "easy_cui.h"
#include "shortest_path/dynamic_pruned_landmark_labeling.h"

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);
  dynamic_pruned_landmark_labeling<64> dpll;
  JLOG_PUT_BENCHMARK("result.construct") dpll.construct(g);
  JLOG_PUT("result.average_label_size",
           (dpll.total_label_num() / g.num_vertices()));
  JLOG_PUT_BENCHMARK("result.query_time_1M") {
    for (size_t i = 0; i < 1000 * 1000; ++i) {
      dpll.query_distance(g, agl::random(g.num_vertices()),
                          agl::random(g.num_vertices()));
    }
  }
  JLOG_PUT_BENCHMARK("result.add_time_1K") {
    for (size_t i = 0; i < 1000; ++i) {
      dpll.add_edge(g, agl::random(g.num_vertices()),
                    agl::random(g.num_vertices()));
    }
  }

  return 0;
}
