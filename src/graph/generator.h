#pragma once
#include "graph/graph.h"
#include "graph/weight_type.h"

namespace agl {
unweighted_edge_list generate_path(V num_vertices);
unweighted_edge_list generate_erdos_renyi(V num_vertices, double avg_deg);
unweighted_edge_list generate_grid(size_t num_rows, size_t num_cols);
unweighted_edge_list generate_barbell(V size_clique);
unweighted_edge_list generate_random_planar(V num_verticese, size_t num_edges);
unweighted_edge_list generate_cycle(V num_vertices);

unweighted_edge_list generate_random_spanning_tree(V num_vertices);

unweighted_edge_list make_undirected(const unweighted_edge_list &es);

template<typename GraphType>
typename GraphType::edge_list_type add_random_weight(const unweighted_edge_list &es);

template<>
unweighted_edge_list add_random_weight<unweighted_graph>(const unweighted_edge_list &es);

template<typename GraphType>
typename GraphType::edge_list_type add_random_weight(const unweighted_edge_list &es) {
  typename GraphType::edge_list_type wes(es.size());
  for (size_t i = 0; i < es.size(); ++i) {
    wes[i].first = es[i].first;
    wes[i].second = typename GraphType::E{es[i].second, random_weight<typename GraphType::W>()};
  }
  return wes;
}

template<typename GraphType>
typename GraphType::edge_list_type add_unit_weight(const unweighted_edge_list &es);

template<>
unweighted_edge_list add_unit_weight<unweighted_graph>(const unweighted_edge_list &es);

template<typename GraphType>
typename GraphType::edge_list_type add_unit_weight(const unweighted_edge_list &es) {
  typename GraphType::edge_list_type wes(es.size());
  for (size_t i = 0; i < es.size(); ++i) {
    wes[i].first = es[i].first;
    wes[i].second = typename GraphType::E{es[i].second, unit_weight<typename GraphType::W>()};
  }
  return wes;
}
}  // namespace agl
