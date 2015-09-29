#include "easy_cui.h"
#include "distance_sketch/distance_sketch.h"
using namespace distance_sketch;

DEFINE_string(ks, "1,8,64", "");
DEFINE_int32(num_retrieval_trials, 10000, "");

namespace {
template<typename Lambda>
void eval(const G &g, size_t k, Lambda lmd) {
  JLOG_PUT("k", k);

  using graph_sketches_type = decltype(lmd(g, k, {}, kFwd));
  unique_ptr<graph_sketches_type> s;

  JLOG_PUT_BENCHMARK("time_construction") {
    s.reset(new graph_sketches_type(lmd(g, k, {}, kFwd)));
  }
  JLOG_PUT("avg_sketch_size", s->average_size());

  {
    double t = get_current_time_sec();
    for (int trial = 0; trial < FLAGS_num_retrieval_trials; ++trial) {
      s->retrieve_sketch(agl::random(g.num_vertices()));
    }
    JLOG_PUT("time_avg_retrieval", (get_current_time_sec() - t) / FLAGS_num_retrieval_trials);
  }
}
}  // namespace

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);

  for (auto k : parse_comma_separated_string<int>(FLAGS_ks)) {
    JLOG_ADD_OPEN("ks") eval(g, k, compute_all_distances_sketches);
  }
  for (auto k : parse_comma_separated_string<int>(FLAGS_ks)) {
    JLOG_ADD_OPEN("ks") eval(g, k, compute_sketch_retrieval_shortcuts_via_ads_unweighted);
  }
  for (auto k : parse_comma_separated_string<int>(FLAGS_ks)) {
    JLOG_ADD_OPEN("ks") eval(g, k, compute_sketch_retrieval_shortcuts_via_ads_fast);
  }

  return 0;
}
