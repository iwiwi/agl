#include "shortest_path/dynamic_pruned_landmark_labeling.h"
using namespace std;
using namespace agl;

int main() {
  unweighted_edge_list es = generate_grid(3, 3);
  G g(es);

  {
    dynamic_pruned_landmark_labeling<0> dpll;
    dpll.construct(g);

    for (V i = 0; i < 9; ++i) {
      for (V j = 0; j < 9; ++j) {
        cout << i << "->" << j << endl;
        cout << dpll.query_distance(g, i, j) << endl;
      }
    }
  }

  return 0;
}