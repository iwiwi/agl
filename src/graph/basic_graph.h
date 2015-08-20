#pragma once
#include "base/base.h"
#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>

namespace agl {
// We don't use scoped enum to avoid writing |static_cast| every time
enum D : uint8_t {
  kFwd = 0,
  kBwd = 1,
  kNumDirections = 2,
};

using direction_underlying_type = std::underlying_type<D>::type;

inline irange<direction_underlying_type> directions() {
  return irange<direction_underlying_type>(kNumDirections);
}

template<typename EdgeType>
class neighbor_range {
 public:
  struct iterator_type {
    typename std::vector<EdgeType>::iterator ite;
    V operator*() const { return to(*ite); }
    void operator++() { ++ite; }
  };

  neighbor_range(std::vector<EdgeType> &edges) :
    i_{edges.begin()}, n_{edges.end()} {}

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
  using W = decltype(weight(E())); // typename E::weight_type;
  using edge_list_type = std::vector<std::pair<V, E>>;

  basic_graph() = default;

  explicit basic_graph(const edge_list_type &es) {
    assign(es);
  }

  void clear();

  void assign(const edge_list_type &es);

  //
  // Graph access
  //
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
    return edges_from_[d][v][i];
  }

  inline size_t degree(V v, D d = kFwd) const {
    return edges_from_[d][v].size();
  }

  inline V num_vertices() const {
    assert(edges_from_[kFwd].size() == edges_from_[kBwd].size());
    return edges_from_[kFwd].size();
  }

  //
  // Dynamic graph update
  //
  template<typename... Args>
  void add_edge(V from, Args&&... args) {
    edges_from_[kFwd][from].emplace_back(std::forward<Args>(args)...);
    const E &e = edges_from_[kFwd][from].back();
    edges_from_[kBwd][to(e)].emplace_back(reverse_edge(from, e));
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
}
}  // namespace agl
