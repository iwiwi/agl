#include "easy_cui.h"
#include "distance_sketch/distance_sketch.h"
using namespace agl::dynamic_graph_update;

int main(int argc, char **argv) {
  G g(easy_cui_init(argc, argv));

  update_scenario<G> s;
  s = generate_scenario_random_addition_and_removal(g.edge_list(), 100);

  {
    distance_sketch::dynamic_all_distances_sketches ads;
    evaluate_scenario(&ads, s);
  }

  {
    distance_sketch::dynamic_sketch_retrieval_shortcuts srs;
    evaluate_scenario(&srs, s);
  }

  return 0;
}
