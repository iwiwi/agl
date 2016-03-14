#pragma once
#include "graph/graph.h"
#include "graph/graph_index_interface.h"

namespace agl {
class dynamic_landmark_distance_estimator
    : public dynamic_graph_index_interface<G>,
      public distance_query_interface<G> {
 public:
  virtual ~dynamic_landmark_distance_estimator() {}
  virtual void construct(const G &g) override;
  virtual W query_distance(const G &g, V v_from, V v_to) override;
};
}  // namespace agl
