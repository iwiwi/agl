#pragma once

#include <unordered_map>
#include "graph/graph.h"
#include "edge_centrality/edge_centrality.h"

DECLARE_string(graphviz_engine);

namespace agl {
class graphviz {
 public:
  using property_map = std::unordered_map<std::string, std::string>;

  //
  // Input
  //
  graphviz() : graphviz_engine_(FLAGS_graphviz_engine) {}
  explicit graphviz(const G &g);
  void read_dot(const std::string &filename);
  void read_dot(std::istream &is);
  void clear();
  void add_edge(V u, V v) {
    edge_property_.emplace(std::make_pair(u, v), property_map());
  }

  //
  // Output
  //
  void output_dot(std::ostream &ofs) const;
  void output_dot(const std::string &dot_filename) const;

  void generate_image(const std::string &img_filename, const std::string &img_type) const;
  void generate_png(const std::string &png_filename) const;
  void generate_eps(const std::string &eps_filename) const;

  void set_graphviz_option(const std::string &graphviz_option) {
    graphviz_option_ = graphviz_option;
  }

  //
  // Direct property manipulation
  //
  void set_graphviz_engine(const std::string &engine) {
    graphviz_engine_ = engine;
  }

  void set_graph_property(const std::string &property, const std::string &value);

  void set_vertex_property(V v,
                           const std::string &property, const std::string value);
  void set_vertex_property(const std::string &property,
                           const std::unordered_map<V, std::string> &value_map);
  void set_vertex_property(const std::string &property,
                           std::function<std::string(V)> value_func);

  void set_edge_property(V u, V v,
                         const std::string &property, const std::string value);
  void set_edge_property(const std::string &property,
                         const std::unordered_map<std::pair<V, V>, std::string> &value_map);
  void set_edge_property(const std::string &property,
                         std::function<std::string(V, V)> value_func);

  std::unordered_map<V, property_map> &vertex_properties() {
    return vertex_property_;
  }

  std::unordered_map<std::pair<V, V>, property_map> &edge_properties() {
    return edge_property_;
  }

  //
  // High-level property configuration
  //
  void merge_reverse_edges();
  void ignore_direction();
  void ignore_isolated_vertex(const G &g);

 private:
  std::string graphviz_engine_, graphviz_option_;

  property_map graph_property_, vertex_common_property_;
  std::unordered_map<V, property_map> vertex_property_;
  std::unordered_map<std::pair<V, V>, property_map> edge_property_;

  static void output_property_string(std::ostream &os, const property_map&);
};

//
// Easy interfaces
//
void graphviz_draw_graph(const G &g,
                         const char *png_filename = "tmp.png",
                         const char *graphviz_engine = FLAGS_graphviz_engine.c_str());

void graphviz_draw_edge_centrality(const G &g,
                                   const edge_centrality_map &ec,
                                   const char *png_filename = "tmp.png",
                                   const char *graphviz_engine = FLAGS_graphviz_engine.c_str());
}  // namespace agl
