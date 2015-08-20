#pragma once

#include "graph/graph.h"

namespace agl {
void graphviz(const G &g, const char *filename, const char *command = "dot");
}  // namespace agl
