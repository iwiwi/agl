#pragma once

#include "graph/graph.h"

DECLARE_string(graphviz_engine);

namespace agl {
void graphviz(const G &g, const char *filename, const char *command = FLAGS_graphviz_engine.c_str());
}  // namespace agl
