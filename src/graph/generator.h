#pragma once
#include "graph/graph.h"
#include "graph/weight_type.h"

namespace agl {
std::vector<std::pair<V, V>> gen_erdos_renyi(V num_vertices, size_t avg_deg);
std::vector<std::pair<V, V>> gen_grid(size_t num_rows, size_t num_cols);

/*
template<typename GraphT>
typename GraphT::edge_list_type add_weight(const unweighted_graph::edge_list_type &es);
*/

inline unweighted_graph::edge_list_type add_weight(const unweighted_graph::edge_list_type &es) {
  return es;
}

template<typename WeightT>
typename weighted_graph<WeightT>::edge_list_type add_weight(const unweighted_graph::edge_list_type &es) {
  typename weighted_graph<WeightT>::edge_list_type wes(es.size());
  for (size_t i = 0; i < es.size(); ++i) {
    wes[i].first = es[i].first;
    wes[i].second.to = to(es[i].second);
    wes[i].second.weight = random_weight<WeightT>();
  }
  return wes;
}
}  // namespace agl
