#pragma once
#include "graph/graph.h"
#include "graph/graph_index_interface.h"

namespace agl {
class dynamic_betweenness_centrality
  : public dynamic_graph_index_interface<G>,
    public betweenness_centrality_query_interface<G> {
 public:
  virtual ~dynamic_betweenness_centrality() {}
  virtual void construct(const G &g) override;
  virtual double query_betweenness_centrality(const GraphType &g, V v) override;
  virtual void add_edge(const G &g, V v_from, const E &e) override;
  virtual void remove_edge(const G &g, V v_from, V v_to) override;
  virtual void remove_vertices(const G &g, V old_num_vertices) override;
  virtual void add_vertices(const G &g, V old_num_vertices) override;

 private:
};
}  // namespace agl
