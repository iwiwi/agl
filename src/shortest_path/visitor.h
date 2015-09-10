#pragma once
#include "heap/heap.h"

// http://pgkiss.web.fc2.com/cxx/parameterize-code.html

namespace agl {
// TODO: specialize for |unweighted_graph|
template<typename GraphType>
class visitor_by_distance {
 public:
  using W = typename GraphType::W;

  visitor_by_distance(const GraphType &g) : g_(g), h_(g.num_vertices()) {}

  template<typename LambdaType>
  void visit(V v_source, LambdaType lmd, D d = kFwd) {
    h_.decrease(v_source, 0);

    while (!h_.empty()) {
      V v = h_.top_vertex();
      W w = h_.top_weight();
      h_.pop();

      if (!lmd(v, w)) continue;

      for (const auto &e : g_.edges(v, d)) {
        h_.decrease(to(e), w + weight(e));
      }
    }

    h_.clear();
  }

 private:
  const GraphType &g_;
  dijkstra_heap<GraphType> h_;
};

template<typename GraphType, typename LambdaType>
void visit_by_distance(const GraphType &g, V v_source, LambdaType lmd, D d = kFwd) {
  visitor_by_distance<GraphType>(g).visit(v_source, lmd, d);
}
}  // namespace agl
