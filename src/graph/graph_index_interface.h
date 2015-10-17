#pragma once

namespace agl {
template<typename GraphType>
class graph_index_interface {
 public:
  virtual ~graph_index_interface() {}
  virtual void construct(const GraphType &g) = 0;
};

template<typename GraphType>
class dynamic_graph_index_interface : public graph_index_interface<GraphType> {
 public:
  using E = typename GraphType::E;
  virtual ~dynamic_graph_index_interface() {}

  virtual void add_edge(const GraphType &g, V v_from, const E &e) = 0;
  virtual void remove_edge(const GraphType &g, V v_from, V v_to) = 0;
  virtual void add_vertices(const GraphType &g, V old_num_vertices) = 0;
  virtual void remove_vertices(const GraphType &g, V old_num_vertices) = 0;
};
}
