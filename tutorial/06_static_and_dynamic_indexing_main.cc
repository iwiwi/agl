// TODO: write explanation
#include "agl.h"
using namespace std;
using namespace agl;

class nothing_index : public graph_dynamic_index_interface<G> {
 public:
  virtual void construct(const G &g) override {
    cout << "CONSTRUCT" << endl;
    pretty_print(g);
  }

  virtual void add_edge(const G &g, V v_from, const E &e) override {
    cout << "NEW EDGE: " << v_from << " -> " << to(e) << endl;
  }

  virtual void remove_edge(const G &g, V v_from, V v_to) override {
    cout << "REMOVED EDGE: " << v_from << "-> " << v_to << endl;
  }

  virtual void add_vertices(const G &g, V old_num_vertices) override {
    cout << "NUM OF VERTICES: " << old_num_vertices << "->" << g.num_vertices() << endl;
  }
  virtual void remove_vertices(const G&, V) override {}
};

int main(int argc, char **argv) {
  // Example of an owned index
  {
    G g(generate_path(3));

    g.construct_observe_and_own_dynamic_index(new nothing_index());

    g.add_edge(0, 2);
    g.remove_edge(0, 2);
    g.add_vertices(8);

    // The |nothing_index| instance is "owned" by g,
    // so it is automatically deleted.
  }

  // Example of a non-owned index
  {
    G g(generate_path(3));

    nothing_index i;
    g.construct_and_observe_dynamic_index(&i);

    g.add_edge(0, 2);
    g.remove_edge(0, 2);
    g.add_vertices(8);
  }

  return 0;
}
