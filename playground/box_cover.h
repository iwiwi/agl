#include "agl.h"
using namespace std;
using namespace agl;

double coverage(const G &g, const vector<V> &s, W rad);

//
// Naive Functions for Tests
//
vector<map<V, V>> naive_build_sketch(const G &g, const W radius, const int k,
                                     const vector<V> &rank,
                                     const vector<bool> &is_covered);
vector<map<V, V>> build_sketch(const G &g, const W radius, const int k,
                               const vector<V> &rank, const vector<V> &inv,
                               const vector<bool> &is_covered);

//
// Radius-based Methods:
//   Returns the set S of selected center nodes.
//   For any vertex v, there is s \in S s.t. d(v, s) <= radius.
//

//! Song et al. 2007 (Section 3.2)
vector<V> box_cover_memb(const G &g, W radius);

//! Schneider et al. 2012
vector<V> box_cover_burning(const G &g, W radius);

//! Akiba et al. 2015
vector<V> box_cover_sketch(const G &g, W radius, const int k,
                           const int path_num);

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
