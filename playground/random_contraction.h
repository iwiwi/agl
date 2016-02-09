#pragma once

DEFINE_bool(prune_if_degree_eq_1, true, "");
DEFINE_int32(solver_iter, 50, "");

struct degree_weight {
  int degree, weight;
  degree_weight() : degree(0), weight(0) {}
  degree_weight(int degree, int weight) : degree(degree), weight(weight) {}
  bool operator<(const degree_weight& r) const {
    return degree < r.degree;
  }
  degree_weight& operator+=(const degree_weight& r) {
    degree += r.degree;
    weight += r.weight;
    return *this;
  }
  degree_weight operator+(const degree_weight& r) const {
    degree_weight ret = *this;
    return ret += r;
  }
};


class weighted_union_find {
  using weight_t = degree_weight;
public:
  weighted_union_find(int n) : n_(n), uf_(n_), par_(n_), weight_(n_) {
    FOR(i, n_) par_[i] = i;
  }
  int root(int x) {
    int uf_root = uf_.root(x);
    return par_[uf_root];
  }

  void unite(int x, int y) {
    x = root(x), y = root(y);
    if (x == y) return;
    uf_.unite(x, y);
    const int more_heavily_v = weight_[x] < weight_[y] ? y : x;
    const int uf_root = uf_.root(more_heavily_v);
    par_[uf_root] = more_heavily_v;
    weight_[root(x)] = weight_[x] + weight_[y];
  }

  bool is_same(int x, int y) {
    return uf_.root(x) == uf_.root(y);
  }

  void add_weight(int x, weight_t dec) {
    weight_[root(x)] += dec;
  }

  weight_t weight(int x) {
    return weight_[root(x)];
  }
private:
  int n_;
  union_find uf_;
  vector<int> par_;
  vector<weight_t> weight_;
};

class random_contraction {

  void build_depth() {
    depth_.resize(sz(binary_tree_edges_), -1);
    parent_weight_.resize(sz(binary_tree_edges_), make_pair(-2, -2));

    FOR(v, num_vertices_) {
      if (depth_[v] >= 0) continue;

      depth_[v] = 0;
      parent_weight_[v] = make_pair(-1, -1);
      queue<int> q;
      q.push(v);
      while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto& to_weight : binary_tree_edges_[u]) {
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

  void contraction(vector<int>& ancestor, vector<unordered_map<V, int>>& contraction_graph_edges, V u, V v, queue<int>& degree_eq_1) {
    assert(!uf_.is_same(u, v));
    V new_vertex = (V)binary_tree_edges_.size();
    u = uf_.root(u), v = uf_.root(v);
    degree_weight uw = uf_.weight(u), vw = uf_.weight(v);
    uf_.unite(u, v);
    if (u != uf_.root(u)) {
      swap(u, v);
      swap(uw, vw);
    }

    int ua = ancestor[u], va = ancestor[v];
    ancestor[u] = new_vertex;

    // 縮約した結果を元にbinary treeの構築を進める
    binary_tree_edges_.emplace_back(); // new_vertex分の確保
    binary_tree_edges_[new_vertex].emplace_back(ua, uw.weight);
    binary_tree_edges_[ua].emplace_back(new_vertex, uw.weight);
    binary_tree_edges_[new_vertex].emplace_back(va, vw.weight);
    binary_tree_edges_[va].emplace_back(new_vertex, vw.weight);

    //縮約した頂点の辺をまとめる
    auto& uset = contraction_graph_edges[u];
    auto& vset = contraction_graph_edges[v];
    int dec_weight = 0;
    for (const auto& edge : vset) {
      V to; int w; tie(to, w) = edge;
      if (to == u) {
        auto it = uset.find(v);
        assert(it != uset.end());
        assert(it->second == w);
        uset.erase(it);
        dec_weight += w * 2;
      } else {
        uset[to] += w;
        auto& toset = contraction_graph_edges[to];
        auto it = toset.find(v);
        toset.erase(it);
        toset[u] += w;
        if (FLAGS_prune_if_degree_eq_1 && sz(toset) == 1) {
          degree_eq_1.push(to);
        }
      }
    }
    vset.clear();
    uf_.add_weight(u, degree_weight(-2, -dec_weight));
  }

public:
  // g の辺を使い、グラフを作成する(勝手にundirectedとして読み替えている)
  random_contraction(G& g) :
    num_vertices_(g.num_vertices()),
    uf_(g.num_vertices()),
    binary_tree_edges_(g.num_vertices()) {
    //unordered -> weight = 1なので、shuffleでよい
    //重みがあるならBITで。(stochastic acceptanceはupdateありだと使えない)
    vector<unordered_map<V, int>> contraction_graph_edges(g.num_vertices());
    typename G::edge_list_type initial_edges(g.edge_list());
    std::shuffle(initial_edges.begin(), initial_edges.end(), agl::random);
    vector<int> ancestor(g.num_vertices());

    FOR(i, num_vertices_) ancestor[i] = i;

    for (auto edge : initial_edges) {
      V u = edge.first, v = to(edge.second);
      contraction_graph_edges[u].emplace(v, 1);
      contraction_graph_edges[v].emplace(u, 1);
    }
    FOR(v, num_vertices_) uf_.add_weight(v, degree_weight(sz(contraction_graph_edges[v]), sz(contraction_graph_edges[v])));

    queue<int> degree_eq_1;
    int pruned_degree_eq_1 = 0;
    if (FLAGS_prune_if_degree_eq_1) {
      FOR(v, num_vertices_) {
        if (sz(contraction_graph_edges[v]) == 1) degree_eq_1.push(v);
      }
    }

    for (const auto& edge : initial_edges) {
      if (FLAGS_prune_if_degree_eq_1) {
        // 次数1の頂点があればそれを縮約する
        while (!degree_eq_1.empty()) {
          V v = degree_eq_1.front(); degree_eq_1.pop();
          if (contraction_graph_edges[v].size() == 0) continue;
          V u;
          tie(u, std::ignore) = *contraction_graph_edges[v].begin();
          contraction(ancestor, contraction_graph_edges, u, v, degree_eq_1);
          if (contraction_graph_edges[u].size() == 1) degree_eq_1.push(u);
          
          pruned_degree_eq_1++;
        }
      }

      V u = edge.first, v = to(edge.second);
      if (uf_.is_same(u, v)) continue;
      contraction(ancestor, contraction_graph_edges, u, v, degree_eq_1);
    }
    // gが連結とは限らないので、uv_costs.size() == g.num_vertices() - 1ではない
    assert(int(binary_tree_edges_.size() - g.num_vertices()) <= g.num_vertices() - 1);

    build_depth();
  }

  int query(V u, V v) const {
    CHECK(u != v);
    CHECK(u < num_vertices_ && v < num_vertices_);
    int ans = numeric_limits<int>::max();
    while (u != v) {
      if (u == -1 || v == -1) return 0; // disconnect
      if (depth_[u] > depth_[v]) {
        ans = min(ans, parent_weight_[u].second);
        u = parent_weight_[u].first;
      } else {
        ans = min(ans, parent_weight_[v].second);
        v = parent_weight_[v].first;
      }
    }
    return ans;
  }

  void single_source_mincut_dfs(V v, V par, int mn, vector<int>& ans) const {
    if (v < num_vertices_) {
      ans[v] = mn;
    }
    for (auto& to : binary_tree_edges_[v]) {
      if (to.first == par) continue;
      int nmn = min(mn, to.second);
      single_source_mincut_dfs(to.first, v, nmn, ans);
    }
  }


  vector<int> single_source_mincut(V u) const {
    vector<int> ans(num_vertices_, numeric_limits<int>::max());
    single_source_mincut_dfs(u, -1, numeric_limits<int>::max(), ans);
    return ans;
  }

private:
  const int num_vertices_;
  weighted_union_find uf_;

  vector<vector<pair<V, int>>> binary_tree_edges_;
  vector<pair<V, int>> parent_weight_;
  vector<int> depth_;
};

class min_cut_query {
public:
  min_cut_query(G& g) {
    fprintf(stderr, "initialize mincut-query tree ...\n");
    for (int i = 0; i < FLAGS_solver_iter; i++) {
      if (i % 10 == 0) fprintf(stderr, "trees .. %d/%d\n", i, FLAGS_solver_iter);
      solvers_.emplace_back(g);
    }
    fprintf(stderr, "completed.\n");
  }

  int query(int u, int v) const {
    int ans = numeric_limits<int>::max();
    for (const auto& solver : solvers_) {
      ans = min(ans, solver.query(u, v));
    }
    return ans;
  }

  vector<int> single_source_mincut(V v) const {
    vector<int> ret;
    for (const auto& solver : solvers_) {
      auto cur_ans = solver.single_source_mincut(v);
      if (sz(ret) == 0) ret = move(cur_ans);
      else {
        FOR(i, sz(ret)) ret[i] = min(ret[i], cur_ans[i]);
      }
    }
    return ret;
  }
private:
  vector<random_contraction> solvers_;
};
