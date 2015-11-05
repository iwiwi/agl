#include "pruned_landmark_labeling.h"
#include "shortest_path/dynamic_pruned_landmark_labeling.h"
using namespace std;
using namespace agl;

int main() {
  V N = 5, M = 3;
  unweighted_edge_list es = generate_ba(N, M);

  size_t E0 = M * (M - 1) / 2;
  size_t E1 = es.size();

  unweighted_edge_list initial_es(es.begin(), es.begin() + E0);
  // Query
  size_t Q = 10;
  vector<pair<V, V>> query(Q);
  for (int i = 0; i < Q; ++i) {
    query[i] = make_pair(agl::random() % E0, agl::random() % E0);
  }

  dynamic_pruned_landmark_labeling<> dpll;
  G g(initial_es);
  dpll.construct(g);

  // UPDATE
  for (int e = E0; e < E1; ++e) {
    dpll.add_edge(g, es[e].first, es[e].second);
  }

  return 0;
}