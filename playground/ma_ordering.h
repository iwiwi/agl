#pragma once

namespace {
class ma_ordering {
public:
  // 入力されるグラフは連結である必要がある
  ma_ordering(const vector<pair<V, V>>& edges, int num_vs)
    : n_(num_vs), root_(0), edges_(n_) {
    for (auto& e : edges) {
      edges_[e.first].push_back(e.second);
      edges_[e.second].push_back(e.first);
    }
  }

  vector<int> run() {
    vector<int> ordering;

    priority_queue<pair<int, V>> pq;
    vector<int> costs(n_);

    auto use_vertex = [&](const V v) {
      if (costs[v] == -1) return;
      costs[v] = -1;
      ordering.push_back(v);
      for (auto to : edges_[v]) {
        if (costs[to] != -1) {
          costs[to]++;
          pq.push(make_pair(costs[to], to));
        }
      }
    };

    use_vertex(root_);
    while (!pq.empty()) {
      auto cost_vertex = pq.top(); pq.pop();
      use_vertex(cost_vertex.second);
    }

    return ordering;
  }

private:
  const int n_;
  V root_;
  vector<vector<V>> edges_;
};
}

vector<V> calc_ma_ordering(const vector<pair<V,V>>& edges, int num_vs) {
  ma_ordering ma(edges, num_vs);
  return ma.run();
}