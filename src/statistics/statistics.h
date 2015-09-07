#pragma once
#include "agl.h"

namespace agl {
std::pair<size_t, size_t> num_triangles_and_wedges(const G &g);
std::vector<std::pair<size_t, size_t>> num_local_triangles_and_wedges(const G &g);

// https://en.wikipedia.org/wiki/Clustering_coefficient#Local_clustering_coefficient
std::vector<double> local_clustering_coefficient(const G &g);
// https://en.wikipedia.org/wiki/Clustering_coefficient#Network_average_clustering_coefficient
double average_clustering_coefficient(const G &g);
// https://en.wikipedia.org/wiki/Clustering_coefficient#Global_clustering_coefficient
double global_clustering_coefficient(const G &g);
}  // namespace agl
