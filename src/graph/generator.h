#pragma once
#include "graph/graph.h"
#include "graph/weight_type.h"

namespace agl {
unweighted_edge_list gen_erdos_renyi(V num_vertices, size_t avg_deg);
unweighted_edge_list gen_grid(size_t num_rows, size_t num_cols);

unweighted_edge_list gen_barbell(V size_clique);
unweighted_edge_list force_undirected(const unweighted_edge_list &es);

/*
template<typename GraphT>
typename GraphT::edge_list_type add_weight(const unweighted_graph::edge_list_type &es);
*/

inline unweighted_edge_list add_weight(const unweighted_edge_list &es) {
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
