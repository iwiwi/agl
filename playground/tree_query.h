#pragma once

class tree_query {

  void build() {
    depth_.resize(num_vertices_, -1);
    parent_weight_.resize(num_vertices_, make_pair(-2, -2));

    FOR(v, num_vertices_) {
      if (depth_[v] >= 0) continue;

      depth_[v] = 0;
      parent_weight_[v] = make_pair(-1, -1);
      queue<int> q;
      q.push(v);
      while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto& to_weight : edges_[u]) {
          int to, weight; tie(to, weight) = to_weight;
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

  void single_source_mincut_dfs(V v, V par, int mn, vector<int>& ans) const {
    ans[v] = mn;
    for (auto& to : edges_[v]) {
      if (to.first == par) continue;
      int nmn = min(mn, to.second);
      single_source_mincut_dfs(to.first, v, nmn, ans);
    }
  }

public:
  static tree_query from_file(const string& path) {
    ifstream ifs(path.c_str());
    return from_file(ifs);
  }

  static tree_query from_file(istream& is) {
    vector<tuple<V, V, int>> input;
    int s, v, weight;
    while (is >> s >> v >> weight) input.emplace_back(s, v, weight);
    return tree_query(std::move(input));
  }

  tree_query(vector<tuple<V, V, int>> input) {
    num_vertices_ = sz(input) + 1;

    edges_.resize(num_vertices_);
    vector<int> depth_;
    for (auto& e : input) {
      int s, v, weight; tie(s, v, weight) = e;
      edges_[s].emplace_back(v, weight);
      edges_[v].emplace_back(s, weight);
    }
    input.clear(); input.shrink_to_fit();

    build();
  }

  vector<int> single_source_mincut(V u) const {
    vector<int> ans(num_vertices_, numeric_limits<int>::max());
    single_source_mincut_dfs(u, -1, numeric_limits<int>::max(), ans);
    return ans;
  }

  int query(V u, V v) const {
    CHECK(u != v);
    CHECK(u < num_vertices_ && v < num_vertices_);
    int ans = numeric_limits<int>::max();
    while (u != v) {
      if (depth_[u] > depth_[v]) swap(u, v);
      ans = min(ans, parent_weight_[v].second);
      v = parent_weight_[v].first;
    }
    return ans;
  }
  
private:
  int num_vertices_;
  vector<pair<V, int>> parent_weight_;
  vector<vector<pair<V, int>>> edges_;
  vector<int> depth_;
};
