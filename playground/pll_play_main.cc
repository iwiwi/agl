#include "shortest_path/dynamic_pruned_landmark_labeling.h"
using namespace std;
using namespace agl;

int main() {
  unweighted_edge_list es = generate_grid(2, 2);
  G g(es);

  {
    dynamic_pruned_landmark_labeling<0> dpll;
    dpll.construct(g);

    for (V i = 0; i < 4; ++i) {
      for (V j = 0; j < 4; ++j) {
        cout << i << "->" << j << " " << dpll.query_distance(g, i, j) << endl;
      }
    }
    cout << dpll.query_distance(g, 0, 3) << endl;
    dpll.add_edge(g, 3, 0);
    cout << dpll.query_distance(g, 0, 3) << endl;

    for (V i = 0; i < 4; ++i) {
      for (V j = 0; j < 4; ++j) {
        cout << i << "->" << j << " " << dpll.query_distance(g, i, j) << endl;
      }
    }
  }

  return 0;
}