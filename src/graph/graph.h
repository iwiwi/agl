#pragma once
#include <cstdint>

namespace agl {
using V = int32_t;
}

#include "direction.h"
#include "edge_type.h"
#include "basic_graph.h"
#include "unweighted_graph.h"
#include "weighted_graph.h"

namespace agl {
using G = unweighted_graph;
using E = unweighted_graph::E;
using W = G::W;
static constexpr W kInfW = std::numeric_limits<W>::max();
}  // namespace agl

#include "basic_graph_io.h"
#include "weight_type.h"
#include "generator.h"
#include "built_in_data.h"
