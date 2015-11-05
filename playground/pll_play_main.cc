#include "shortest_path/dynamic_pruned_landmark_labeling.h"
using namespace std;
using namespace agl;

int main() {
  unweighted_edge_list es = generate_grid(3, 3);
  G g(es);

  {
    dynamic_pruned_landmark_labeling<> dpll;
    dpll.construct(g);
    cout << dpll.query_distance(g, 0, 5) << endl;
    dpll.add_edge(g, 0, 5);
    cout << dpll.query_distance(g, 0, 5) << endl;
  }

  return 0;
}