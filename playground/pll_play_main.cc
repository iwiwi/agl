#include "pruned_landmark_labeling.h"
using namespace std;
using namespace agl;

int main() {
  unweighted_edge_list es = generate_grid(3, 3);
  G g(es);

  {
    PrunedLandmarkLabeling<0> pll;
    pll.ConstructIndex(es);
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < 9; ++j) {
        cout << i << "->" << j << endl;
        cout << pll.QueryDistance(i, j) << endl;
      }
    }
  }

  return 0;
}