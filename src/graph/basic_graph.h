#pragma once
#include "base/base.h"
#include "direction.h"
#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <algorithm>

namespace agl {
template<typename EdgeType>
class neighbor_range {
 public:
  struct iterator_type {
    typename std::vector<EdgeType>::const_iterator ite;
    V operator*() const { return to(*ite); }
    void operator++() { ++ite; }
    bool operator!=(const iterator_type &i) const { return ite < i.ite; }
  };

  neighbor_range(const std::vector<EdgeType> &edges) :
    i_{edges.begin()}, n_{edges.end()} {}

  iterator_type begin() { return i_; }
  iterator_type end() { return n_; }

 private:
  iterator_type i_, n_;
};

template<typename EdgeType>
class undirected_neighbor_range {
 public:
  struct iterator_type {
    typename std::vector<EdgeType>::const_iterator i0, i1;
    const undirected_neighbor_range &p;

    V operator*() const {
      if (i0 == p.n_.i0) return to(*i1);
      if (i1 == p.n_.i1) return to(*i0);
      return std::min(to(*i0), to(*i1));
    }

    bool operator!=(const iterator_type &i) const {
      return std::tie(i0, i1) != std::tie(i.i0, i.i1);
    }

    void operator++() {
      if (i0 == p.n_.i0) ++i1;
      else if (i1 == p.n_.i1) ++i0;
      else if (to(*i0) == to(*i1)) {
        ++i0;
        ++i1;
      } else {
        ++(to(*i0) < to(*i1) ? i0 : i1);
      }
    }
  };

  undirected_neighbor_range(const std::vector<EdgeType> &es0,
                            const std::vector<EdgeType> &es1)
  : i_{es0.begin(), es1.begin(), *this}, n_{es0.end(), es1.end(), *this} {}

  iterator_type begin() { return i_; }
  iterator_type end() { return n_; }

 private:
  iterator_type i_, n_;
};

template<typename EdgeType>
class basic_graph {
public:
  using V = agl::V;
  using E = EdgeType;
  using W = decltype(weight(E()));
  using edge_list_type = std::vector<std::pair<V, E>>;

  basic_graph() = default;

  explicit basic_graph(const edge_list_type &es) {
    assign(es);
  }

  void clear();

  void assign(const edge_list_type &es);

  edge_list_type edge_list(D d = kFwd) const;

  //
  // Graph access
  //
  inline V num_vertices() const {
    assert(edges_from_[kFwd].size() == edges_from_[kBwd].size());
    return edges_from_[kFwd].size();
  }

  inline size_t num_edges() const {
    size_t n = 0;
    for (V v : vertices()) n += degree(v);
    return n;
  }

  inline irange<V> vertices() const {
    return make_irange(num_vertices());
  }

  inline const std::vector<E> &edges(V v, D d = kFwd) const {
    return edges_from_[d][v];
  }

  inline const E &edge(V v, size_t i, D d = kFwd) const {
    return edges_from_[d][v][i];
  }

  inline neighbor_range<E> neighbors(V v, D d = kFwd) const {
    return neighbor_range<E>(edges_from_[d][v]);
  }

  inline V neighbor(V v, size_t i, D d = kFwd) const {
    return to(edges_from_[d][v][i]);
  }

  inline size_t degree(V v, D d = kFwd) const {
    return edges_from_[d][v].size();
  }

  //
  // Dynamic graph update
  //
  template<typename... Args>
  void add_edge(V v_from, Args&&... args) {
    edges_from_[kFwd][v_from].emplace_back(std::forward<Args>(args)...);
    const E &e = edges_from_[kFwd][v_from].back();
    const V v_to = to(e);
    edges_from_[kBwd][v_to].emplace_back(reverse_edge(v_from, e));

    // TODO: faster by insertion
    auto cmp = [](const E &e0, const E &e1) { return to(e0) < to(e1); };
    std::sort(edges_from_[kFwd][v_from].begin(), edges_from_[kFwd][v_from].end(), cmp);
    std::sort(edges_from_[kBwd][v_to].begin(), edges_from_[kBwd][v_to].end(), cmp);
  }

private:
  std::vector<std::vector<E>> edges_from_[kNumDirections];
};

namespace internal {
template<typename EdgePairType>
V num_vertices(const std::vector<EdgePairType> &es) {
  V num_vs = 0;
  for (const auto &p : es) {
    num_vs = std::max(num_vs, std::max(p.first, to(p.second)) + 1);
  }
  return num_vs;
}
}  // namespace internal

//
// Member functions of |basic_graph|
//
template<typename EdgeType>
void basic_graph<EdgeType>::assign(const basic_graph<EdgeType>::edge_list_type &es) {
  V num_vs = internal::num_vertices(es);
  edges_from_[kFwd].assign(num_vs, {});
  edges_from_[kBwd].assign(num_vs, {});
  for (const auto &p : es) {
    edges_from_[kFwd][p.first].emplace_back(p.second);
    edges_from_[kBwd][to(p.second)].emplace_back(reverse_edge(p.first, p.second));
  }
  for (D d : directions()) {
    for (V v : vertices()) {
      auto cmp = [](const EdgeType &e0, const EdgeType &e1) {
        return to(e0) < to(e1);
      };
      std::sort(edges_from_[d][v].begin(), edges_from_[d][v].end(), cmp);
    }
  }
}

template<typename EdgeType>
typename basic_graph<EdgeType>::edge_list_type basic_graph<EdgeType>::edge_list(D d) const {
  decltype(edge_list(d)) el;
  for (V v : vertices()) {
    for (const E &e : edges(v, d)) {
      el.emplace_back(v, e);
    }
  }
  return el;
}

//
// Easy utility functions
//
template<typename GraphType>
bool is_adjacent(const GraphType &g, V u, V v, D d = kFwd) {
  const auto &r = g.edges(u, d);
  return std::binary_search(r.begin(), r.end(), v);
}

template<typename GraphType>
inline undirected_neighbor_range<typename GraphType::E> undirected_neighbors(const GraphType &g, V v) {
  return undirected_neighbor_range<typename GraphType::E>(g.edges(v, kFwd), g.edges(v, kBwd));
}
}  // namespace agl
