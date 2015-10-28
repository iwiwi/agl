#include "easy_cui.h"
#include "distance_sketch/distance_sketch.h"

int main(int argc, char **argv) {
  dynamic_index_evaluation_scenario<G> s;
  s.initial_graph = easy_cui_init(argc, argv);
  s.add_workload_edge_addition_and_removal_random(100);

  distance_sketch::dynamic_all_distances_sketches ads;
  s.evaluate((dynamic_graph_index_interface<G>*)&ads);
  return 0;  // Only ADS

  distance_sketch::dynamic_sketch_retrieval_shortcuts srs;
  s.evaluate((dynamic_graph_index_interface<G>*)&srs);
}
