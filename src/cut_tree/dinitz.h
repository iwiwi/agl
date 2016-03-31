#pragma once
#include <base/base.h>
#include <graph/graph.h>

namespace agl {
namespace cut_tree_internal {
class dinitz {
  struct E {
    int to_, rev_, cap_;
    E(int to, int rev_, int cap_) : to_(to), rev_(rev_), cap_(cap_) {}
  };

  void bfs(int s);
  int dfs(int v, int t, int f);
  void add_undirected_edge(int f, int t, int c);
  void reset_graph();

public:
  dinitz(const G& g);
  int max_flow(int s, int t);

  std::vector<E>& edges(V v) { return e_[v]; }
  V to(const E& e) const { return e.to_; }
  int cap(E& e) { return e.cap_; }
  E& rev(const E& e_in) { return e_[e_in.to_][e_in.rev_]; }

private:
  G g;
  std::vector<int> level_, iter_;
  std::vector<std::vector<E>> e_;
};

} //namespace cut_tree_internal
} //namespace agl