// This is how we evaluate dynamic indexing methods
// TODO: write explanation
#include "agl.h"
using namespace std;
using namespace agl;

// An example of dynamic indices
class my_slow_index : public graph_dynamic_index_interface<G> {
 public:
  virtual void construct(const G &g) override {
    cout << "CONSTRUCT" << endl;
    usleep(300000);
  }

  virtual void add_edge(const G &g, V v_from, const E &e) override {
    cout << "NEW EDGE: " << v_from << " -> " << to(e) << endl;
    usleep(100000);
  }

  virtual void remove_edge(const G &g, V v_from, V v_to) override {
    cout << "REMOVED EDGE: " << v_from << "-> " << v_to << endl;
    usleep(200000);
  }

  virtual void add_vertices(const G &g, V old_num_vertices) override {}
  virtual void remove_vertices(const G&, V) override {}
};

int main() {
  // Prepare the scenario
  dynamic_index_evaluation_scenario<G> s;
  s.initial_graph.assign(generate_path(10));
  s.add_workload_edge_addition_and_removal_random(10);

  // Evaluate
  my_slow_index i;
  s.evaluate(&i);
}
