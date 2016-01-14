#include "agl.h"
using namespace std;
using namespace agl;

class coverage_manager {
 private:
  V num_v;
  vector<W> dist;
  V cnt = 0;
  W radius;
  double goal_coverage;
  coverage_manager() {}

 public:
  coverage_manager(const G &g, const W r, const double c) {
    num_v = g.num_vertices();
    dist.assign(num_v, num_v);
    radius = r;
    assert(c <= 1.0);
    goal_coverage = c;
  }

  void add(const G &g, V new_v) {
    queue<V> que;
    if (dist[new_v] == num_v) cnt++;
    dist[new_v] = 0;
    que.push(new_v);

    for (W d = 0; d < radius; ++d) {
      size_t s = que.size();
      for (size_t t = 0; t < s; t++) {
        V a = que.front();
        que.pop();
        for (V neighbor : g.neighbors(a)) {
          if (dist[neighbor] <= d + 1) continue;
          if (dist[neighbor] == num_v) cnt++;
          dist[neighbor] = d + 1;
          que.push(neighbor);
        }
      }
    }
  }

  double get_current_coverage() { return (double)cnt / num_v; }

  bool is_covered() { return get_current_coverage() >= goal_coverage; }
};

vector<pair<W, V>> find_analytical_solution(const string &type, V u, V v,
                                            const G &g);

double naive_coverage(const G &g, const vector<V> &s, W rad);
double coverage(const G &g, const vector<V> &s, W rad);

//
// Naive Functions for Tests
//
vector<vector<V>> naive_build_sketch(const G &g, const W radius, const int k,
                                     const vector<V> &rank,
                                     const vector<V> &inv,
                                     const vector<bool> &is_covered);
vector<vector<V>> build_sketch(const G &g, const W radius, const int k,
                               const vector<V> &rank, const vector<V> &inv,
                               const vector<bool> &is_covered);
vector<vector<V>> build_sketch(const G &g, const W radius, const int k,
                               const vector<V> &rank, const vector<V> &inv,
                               const vector<bool> &is_covered, bool &use_memb,
                               size_t size_upper_bound);

void naive_select_greedily(const G &g, const vector<vector<V>> &X,
                           vector<V> &centers, vector<bool> &centered,
                           const int k);
void select_greedily(const G &g, const vector<vector<V>> &X, vector<V> &centers,
                     vector<bool> &centered, const int k, coverage_manager &cm);
void select_lazy_greedily(const G &g, const vector<vector<V>> &X,
                          vector<V> &centers, vector<bool> &centered,
                          coverage_manager &cm);

//
// Radius-based Methods:
//   Returns the set S of selected center nodes.
//   For any vertex v, there is s \in S s.t. d(v, s) <= radius.
//

//! Song et al. 2007 (Section 3.2)
vector<V> box_cover_memb(const G &g, W radius);

//! Schneider et al. 2012
vector<V> box_cover_burning(const G &g, W radius);

//! Akiba et al. 2016
vector<V> box_cover_sketch(const G &g, W radius, const int k,
                           const int pass_num, double &aim_coverage,
                           double lazy = 1.0);

//
// Diameter-based Methods:
//   returns the sets of vertices with the limited diameter.
//   Any vertex is covered by a set.
//

//! Song et al. 2005
vector<vector<V>> box_cover_original(const G &g, W diameter);

//! Song et al. 2007 (Section 3.1)
vector<vector<V>> box_cover_cbb(const G &g, W diameter);

//! Song et al. 2007 (Section 2)
vector<vector<V>> box_cover_coloring(const G &g, W diameter);
