#include "shortest_path/dynamic_pruned_landmark_labeling.h"
using namespace std;
using namespace agl;

int main() {
  int N = 3;
  unweighted_edge_list es = generate_grid(N, N);
  G g(es);

  {
    dynamic_pruned_landmark_labeling<0> dpll;
    dpll.construct(g);

    cout << dpll.query_distance(g, 0, 4) << endl;
    dpll.add_edge(g, 0, 4);
    cout << dpll.query_distance(g, 0, 4) << endl;

    
  }

  return 0;
}