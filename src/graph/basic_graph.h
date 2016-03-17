#pragma once
#include "base/base.h"
#include "direction.h"
#include "graph_index_interface.h"
#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <memory>

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
  using G = basic_graph<EdgeType>;
  using V = agl::V;
  using E = EdgeType;
  using W = decltype(weight(E()));
  using edge_list_type = std::vector<std::pair<V, E>>;

  basic_graph() = default;

  explicit basic_graph(const edge_list_type &es, V num_vs = -1) {
    assign(es, num_vs);
  }

  basic_graph(const basic_graph<EdgeType> &g) : edges_from_{g.edges_from_[0], g.edges_from_[1]} {
    CHECK_MSG(g.graph_indices_.empty(), "Graphs with attached indices cannot be copied");
    CHECK_MSG(g.graph_dynamic_indices_.empty(), "Graphs with attached indices cannot be copied");
  }

  basic_graph<EdgeType> &operator=(const basic_graph<EdgeType> &g) {
    CHECK_MSG(g.graph_indices_.empty(), "Graphs with attached indices cannot be copied");
    CHECK_MSG(g.graph_dynamic_indices_.empty(), "Graphs with attached indices cannot be copied");
    edges_from_[0] = g.edges_from_[0];
    edges_from_[1] = g.edges_from_[1];
    return *this;
  }

  void clear();

  void assign(const edge_list_type &es, V num_vs = -1);
  void assign(std::vector<std::vector<E>> edges_from);

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

    // Notification to dynamic indices
    for (auto i : graph_dynamic_indices_) {
      i->add_edge(*this, v_from, e);
    }
  }

  void remove_edge(V v_from, V v_to) {
    {
      auto &es = edges_from_[kFwd][v_from];
      auto i = std::remove_if(es.begin(), es.end(), [&](const E &e) { return to(e) == v_to; });
      es.erase(i, es.end());
    }
    {
      auto &es = edges_from_[kBwd][v_to];
      auto i = std::remove_if(es.begin(), es.end(), [&](const E &e) { return to(e) == v_from; });
      es.erase(i, es.end());
    }

    // Notification to dynamic indices
    for (auto i : graph_dynamic_indices_) {
      i->remove_edge(*this, v_from, v_to);
    }
  }

  void add_vertices(V new_num_vertices) {
    V old_num_vertices = num_vertices();
    CHECK(new_num_vertices >= old_num_vertices);
    edges_from_[kFwd].resize(new_num_vertices);
    edges_from_[kBwd].resize(new_num_vertices);

    // Notification to dynamic indices
    for (auto i : graph_dynamic_indices_) {
      i->add_vertices(*this, old_num_vertices);
    }
  }

  void remove_vertices(V new_num_vertices) {
    V old_num_vertices = num_vertices();
    CHECK(new_num_vertices <= old_num_vertices);
    for (V v = new_num_vertices; v < old_num_vertices; ++v) {
      CHECK(degree(v, kFwd) == 0);
      CHECK(degree(v, kBwd) == 0);
    }

    edges_from_[kFwd].resize(new_num_vertices);
    edges_from_[kBwd].resize(new_num_vertices);

    // Notification to dynamic indices
    for (auto i : graph_dynamic_indices_) {
      i->remove_vertices(*this, old_num_vertices);
    }
  }

  //
  // Indexing
  //
  template<typename GraphIndexType>
  void construct_and_own_index(GraphIndexType *index) {
    index->construct(*this);  // Construct
    graph_indices_.emplace_back(index);  // Own
  }

  template<typename GraphIndexType>
  void construct_and_observe_dynamic_index(GraphIndexType *index) {
    index->construct(*this);  // Construct
    graph_dynamic_indices_.emplace_back(index);  // Observe
  }

  template<typename GraphIndexType>
  void construct_observe_and_own_dynamic_index(GraphIndexType *index) {
    index->construct(*this);  // Construct
    graph_dynamic_indices_.emplace_back(index);  // Observe
    graph_indices_.emplace_back(index);  // Own
  }

private:
  std::vector<std::vector<E>> edges_from_[kNumDirections];

  std::vector<std::unique_ptr<graph_index_interface<G>>> graph_indices_;
  std::vector<dynamic_graph_index_interface<G>*> graph_dynamic_indices_;  // Observer pattern
};

//
// Member functions of |basic_graph|
//
template<typename EdgePairType>
V num_vertices_from_edge_list(const std::vector<EdgePairType> &es) {
  V num_vs = 0;
  for (const auto &p : es) {
    num_vs = std::max(num_vs, std::max(p.first, to(p.second)) + 1);
  }
  return num_vs;
}

template<typename EdgeType>
void basic_graph<EdgeType>::assign(const typename basic_graph<EdgeType>::edge_list_type &es, V num_vs) {
  if (num_vs == -1) {
    num_vs = num_vertices_from_edge_list(es);
  }
  edges_from_[kFwd].assign(num_vs, {});
  edges_from_[kBwd].assign(num_vs, {});
  for (const auto &p : es) {
    assert(p.first < num_vs);
    assert(to(p.second) < num_vs);
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

  // TODO: reconstruction?
  assert(graph_indices_.empty());
  assert(graph_dynamic_indices_.empty());
}

template<typename EdgeType>
void basic_graph<EdgeType>::assign(std::vector<std::vector<E>> sorted_edges_from) {
  V num_vs = static_cast<V>(sorted_edges_from.size());
  edges_from_[kFwd] = std::move(sorted_edges_from);
  edges_from_[kBwd].assign(num_vs, {});
  for(V v = 0; v < num_vs; v++) {
    for(size_t i = 0; i < edges_from_[kFwd][v].size(); i++) {
      const auto& e = edges_from_[kFwd][v][i];
      if(i != 0) {
        //is sorted?
        const auto& prev_e = edges_from_[kFwd][v][i-1];
        CHECK(to(prev_e) < to(e));
      }
      edges_from_[kBwd][to(e)].emplace_back(reverse_edge(v, e));
    }
  }

  assert(graph_indices_.empty());
  assert(graph_dynamic_indices_.empty());
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
