#include "pruned_landmark_labeling.h"
#include "shortest_path/dynamic_pruned_landmark_labeling.h"
using namespace std;
using namespace agl;

int main() {
  unweighted_edge_list es = generate_grid(3, 3);
  G g(es);
  pretty_print(g);

  dynamic_pruned_landmark_labeling<> dpll;
  dpll.construct(g);
  cout << dpll.query_distance(g, 0, 5);
  // UPDATE
  dpll.add_edge(g, 0, 5);
  cout << dpll.query_distance(g, 0, 5);

  return 0;
}