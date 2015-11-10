#include "shortest_path/dynamic_pruned_landmark_labeling.h"
using namespace std;
using namespace agl;

int main() {
  unweighted_edge_list es = {{0, 1}, {1, 2}};
  G g(es);

  {
    dynamic_pruned_landmark_labeling<0> dpll;
    dpll.construct(g);

    cout << dpll.query_distance(g, 0, 2) << endl;
    dpll.add_edge(g, 0, 2);
    cout << dpll.query_distance(g, 0, 2) << endl;
    // dpll.add_edge(g, 5, 0);
    // cout << dpll.query_distance(g, 5, 0) << endl;
  }

  return 0;
}