#include "agl.h"

namespace box_cover_interface {
class coverage_manager {
 private:
  agl::V num_v;
  std::vector<agl::W> dist;
  agl::V cnt;
  agl::W radius;
  double least_coverage;
  coverage_manager() {}

 public:
  /**
   * Coverage Manager
   * \param g is the graph to cover
   * \param r is the radius of each box
   * \param c is the least coverage
   */
  coverage_manager(const agl::G &g, const agl::W r, const double c)
      : cnt(0), radius(r), least_coverage(c) {
    assert(least_coverage <= 1.0);
    num_v = g.num_vertices();
    dist.assign(num_v, num_v);
  }

  void add(const agl::G &g, const agl::V new_v) {
    std::queue<agl::V> que;
    if (dist[new_v] == num_v) cnt++;
    dist[new_v] = 0;
    que.push(new_v);

    for (agl::W d = 0; d < radius; ++d) {
      size_t s = que.size();
      for (size_t t = 0; t < s; t++) {
        agl::V a = que.front();
        que.pop();
        for (agl::V neighbor : g.neighbors(a)) {
          if (dist[neighbor] <= d + 1) continue;
          if (dist[neighbor] == num_v) cnt++;
          dist[neighbor] = d + 1;
          que.push(neighbor);
        }
      }
    }
  }
  double get_current_coverage() { return (double)cnt / num_v; }
  bool is_covered() { return get_current_coverage() >= least_coverage; }
  bool v_covered(agl::V v) const { return dist[v] <= radius; }
  bool is_center(agl::V v) const { return dist[v] == 0; }
};

double naive_coverage(const agl::G &g, const std::vector<agl::V> &s, agl::W rad);
std::vector<std::pair<agl::W, agl::V>> find_analytical_solution(const std::string &type, agl::V u, agl::V v, const agl::G &g);

//
// Naive Functions for Tests
//
std::vector<std::vector<agl::V>> naive_build_sketch(const agl::G &g, const agl::W radius, const int k, const std::vector<agl::V> &rank, const std::vector<agl::V> &inv, const std::vector<bool> &is_covered);
std::vector<std::vector<agl::V>> build_sketch(const agl::G &g, const agl::W radius, const int k, const std::vector<agl::V> &rank, const std::vector<agl::V> &inv, const coverage_manager &cm);
std::vector<std::vector<agl::V>> build_sketch(const agl::G &g, const agl::W radius, const int k, const std::vector<agl::V> &rank, const std::vector<agl::V> &inv, const coverage_manager &cm, bool &use_memb, size_t index_size_limit);
void naive_select_greedily(const agl::G &g, const std::vector<std::vector<agl::V>> &X, std::vector<agl::V> &centers, std::vector<bool> &centered, const int k);
void select_greedily(const agl::G &g, const std::vector<std::vector<agl::V>> &X, std::vector<agl::V> &centers, const int k, coverage_manager &cm);
void select_lazy_greedily(const agl::G &g, const std::vector<std::vector<agl::V>> &X, const std::vector<agl::V> &rank, const std::vector<agl::V> &inv, std::vector<agl::V> &centers, coverage_manager &cm);
}  // namespace box_cover_interface

namespace agl {
//
// Radius-based Methods:
//   Returns the set S of selected center nodes.
//   For any vertex v, there is s \in S s.t. d(v, s) <= radius.
//

//! Song et al. 2007 (Section 3.2)
std::vector<V> box_cover_memb(const G &g, W radius);

//! Schneider et al. 2012
std::vector<V> box_cover_burning(const G &g, W radius);

//! Akiba et al. 2016
std::vector<V> box_cover_sketch(const G &g, W radius, const int k, const int pass, double least_coverage = 1.0, double alpha = 1.0);
std::vector<V> box_cover_sketch(const G &g, W radius, const int k, const int pass, box_cover_interface::coverage_manager &cm, double alpha = 1.0);

//
// Diameter-based Methods:
//   returns the sets of vertices with the limited diameter.
//   Any vertex is covered by a set.
//

//! Song et al. 2007 (Section 3.1)
std::vector<V> box_cover_cbb(const G &g, W diameter);

//! Song et al. 2007 (Section 2)
std::vector<std::pair<W, size_t>> box_cover_coloring(const G &g, W diameter);
}  // namespace agl