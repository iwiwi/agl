// TODO: write explanation
#include "agl.h"
#include "graph/dynamic_index_evaluation.h"
using namespace std;
using namespace agl;

class my_slow_index : public graph_dynamic_index_interface<G> {
 public:
  virtual void construct(const G &g) override {
    cout << "CONSTRUCT" << endl;
    pretty_print(g);
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
  unweighted_edge_list es = generate_path(10);
  dynamic_update_scenario<G> s;

  s.initial_graph.assign(es);
  dynamic_update_scenario<G>::generator::add_workload_edge_addition_and_removal_random(&s, 10);

  my_slow_index i;
  dynamic_update_scenario<G>::evaluator::evaluate(s, &i);
}
