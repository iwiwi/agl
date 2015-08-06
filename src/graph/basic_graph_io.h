#pragma once
#include "base/base.h"
#include "graph.h"

namespace agl {
template<typename GraphType = G>
GraphType read_graph_tsv(std::istream &is = std::cin) {
  using E = typename GraphType::E;

  V v;
  E e;
  typename GraphType::edge_list_type es;
  while (is >> v >> e) es.emplace_back(v, e);
  return GraphType(es);
}

template<typename GraphType = G>
GraphType read_graph_tsv(const char *filename) {
  if (strcmp(filename, "-") == 0) {
    return read_graph_tsv(std::cin);
  } else {
    std::ifstream ifs(filename);
    CHECK_PERROR(ifs);
    return read_graph_tsv<GraphType>(ifs);
  }
}

template<typename GraphT = G>
void write_graph_tsv(const GraphT &g, std::ostream &os = std::cout) {
  for (auto v : g.vertices()) {
    for (auto e : g.edges(v)) {
      os << e << std::endl;
    }
  }
}

template<typename GraphType = G>
void write_graph_tsv(const GraphType &g, const char *filename) {
  if (strcmp(filename, "-") == 0) {
    return write_graph_tsv(g, std::cout);
  } else {
    std::ifstream ofs(filename);
    CHECK_PERROR(ofs);
    return write_graph_tsv(g, ofs);
  }
}

template<typename GraphType = G>
void pretty_print(const GraphType &g, std::ostream &os = std::cerr, D d = kFwd) {
  for (V v = 0; v < std::min(10, g.num_vertices()); ++v) {
    os << v << " -> ";
    size_t i = 0;
    for (; i < std::min(size_t(10), g.degree(v, d)); ++i) {
      if (i > 0) os << ", ";
      os << g.edge(v, i, d);
    }
    if (i < g.degree(v, d)) os << "...";
    os << std::endl;
  }
}
}  // namespace agl
