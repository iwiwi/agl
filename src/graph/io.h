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
void pretty_print(const GraphType &g, std::ostream &os = std::cerr) {
  static constexpr V kLimitNumVertices = 5;
  static constexpr size_t kLimitNumEdges = 10;

  os << "=========" << std::endl;
  os << "  Vertices: " << g.num_vertices() << std::endl;
  os << "  Edges: " << g.num_edges() << std::endl;
  os << "  Type: " << typename_of(g) << std::endl;
  for (D d : directions()) {
    os << "----------" << std::endl;
    for (V v = 0; v < std::min(kLimitNumVertices, g.num_vertices()); ++v) {
      os << "  " << v << (d == kFwd ? " -> " : " <- ");
      for (size_t i = 0; i < std::min(kLimitNumEdges, g.degree(v, d)); ++i) {
        if (i > 0) os << ", ";
        os << g.edge(v, i, d);
      }
      if (kLimitNumEdges < g.degree(v, d)) os << ", ...";
      os << std::endl;
    }
    if (kLimitNumVertices < g.num_vertices()) {
      os << "  ..." << std::endl;
    }
  }
  os << "=========" << std::endl;
}
}  // namespace agl
