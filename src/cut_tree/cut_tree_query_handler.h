#pragma once
#include <base/base.h>
#include <graph/graph.h>
#include <vector>
#include <queue>
#include <utility>
#include <string>
#include <fstream>

namespace agl {
class cut_tree_query_handler {
  void build(std::vector<std::vector<std::pair<V, int>>>& edges) {
    depth_.resize(num_vertices_, -1);
    parent_weight_.resize(num_vertices_, make_pair(-2, -2));

    for (V v = 0; v < num_vertices_; v++) {
      if (depth_[v] >= 0) continue;

      depth_[v] = 0;
      parent_weight_[v] = std::make_pair(-1, -1);
      std::queue<int> q;
      q.push(v);
      while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto& to_weight : edges[u]) {
          int to, weight; std::tie(to, weight) = to_weight;
          if (depth_[to] >= 0) continue;
          depth_[to] = depth_[u] + 1;
          parent_weight_[to].first = u;
          parent_weight_[to].second = weight;
          // printf("%d - %d\n",u, to);
          q.push(to);
        }
      }
    }
  }

public:
  static cut_tree_query_handler from_file(const std::string& path) {
    std::ifstream ifs(path.c_str());
    CHECK(ifs.is_open());
    return from_file(ifs);
  }

  static cut_tree_query_handler from_file(std::istream& is) {
    std::vector<std::tuple<V, V, int>> input;
    int s, v, weight;
    while (is >> s >> v >> weight) input.emplace_back(s, v, weight);
    return cut_tree_query_handler(std::move(input));
  }
  cut_tree_query_handler() {}

  cut_tree_query_handler(std::vector<std::tuple<V, V, int>> input) {
    num_vertices_ = (int)input.size() + 1;

    std::vector<std::vector<std::pair<V, int>>> edges(num_vertices_);
    std::vector<int> depth_;
    for (auto& e : input) {
      int s, v, weight; std::tie(s, v, weight) = e;
      edges[s].emplace_back(v, weight);
      edges[v].emplace_back(s, weight);
    }
    input.clear(); input.shrink_to_fit();

    build(edges);
  }

  int query(V u, V v) const {
    CHECK(u != v);
    CHECK(u < num_vertices_ && v < num_vertices_);
    int ans = std::numeric_limits<int>::max();
    while (u != v) {
      if (depth_[u] > depth_[v]) std::swap(u, v);
      ans = std::min(ans, parent_weight_[v].second);
      v = parent_weight_[v].first;
    }
    return ans;
  }

  const int num_vertices() { return num_vertices_; }

public:
  int num_vertices_;
  std::vector<std::pair<V, int>> parent_weight_;
  std::vector<int> depth_;
};
} // namespace agl