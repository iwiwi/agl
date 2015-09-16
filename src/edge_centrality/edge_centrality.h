#pragma once
#include <unordered_map>
#include "graph/graph.h"

DECLARE_string(edge_betweenness_centrality_algorithm);
DECLARE_int32(edge_betweenness_centrality_num_samples);

namespace agl {
using edge_centrality_map = std::unordered_map<std::pair<V, V>, double>;
edge_centrality_map merge_edge_centrality_map_entries_for_undirected_graph
(const edge_centrality_map &ecm);

edge_centrality_map edge_betweenness_centrality(const G &g);
edge_centrality_map edge_betweenness_centrality_naive(const G &g);
edge_centrality_map edge_betweenness_centrality_sample_slow(const G &g);
edge_centrality_map edge_betweenness_centrality_sample(const G &g);
}  // namespace agl
