#include "shortest_path/dynamic_pruned_landmark_labeling.h"
using namespace std;
using namespace agl;

int main() {
  unweighted_edge_list es = generate_grid(2, 2);
  G g(es);
  pretty_print(g);

  {
    dynamic_pruned_landmark_labeling<0> dpll;
    dpll.construct(g);

    for (V i = 0; i < 4; ++i) {
      for (V j = 0; j < 4; ++j) {
        cout << i << "->" << j << " " << dpll.query_distance(g, i, j) << endl;
      }
    }
  }

  return 0;
}