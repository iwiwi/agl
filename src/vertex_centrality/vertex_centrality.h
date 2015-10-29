#pragma once
#include "graph/graph.h"

namespace agl {
std::vector<double> vertex_centrality_degree(const G &g);

/// 1 / (sum of distances to all vertices), higher is better, currently O(nm)
std::vector<double> vertex_centrality_classic_closeness(const G &g);
}  // namespace agl
