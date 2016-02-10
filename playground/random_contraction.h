#pragma once
#include "ConnectedComponentsFilter.h"

DEFINE_bool(prune_if_degree_eq_1, true, "");
DEFINE_int32(solver_iter, 50, "");

template<class weight_t>
class weighted_union_find {
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

class tree_minimum_query {
  void build(const vector<tuple<V, V, int>>& edges) {
    vector<vector<pair<V, int>>> g(n_);
    for (auto& e : edges) {
      V u, v; int w; tie(u, v, w) = e;
      g[u].emplace_back(v, w);
      g[v].emplace_back(u, w);
    }

    depth_.resize(n_, -1);
    parent_weight_.resize(n_, make_pair(-2, -2));

    FOR(v, n_) {
      if (depth_[v] >= 0) continue;

      depth_[v] = 0;
      parent_weight_[v] = make_pair(-1, -1);
      queue<int> q;
      q.push(v);
      while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto& to_weight : g[u]) {
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

public:
  tree_minimum_query() = default;
  tree_minimum_query(const vector<tuple<V, V, int>>& edges, int num_vs) : n_(num_vs) {
    build(edges);
  }

  int query(V u, V v) const {
    CHECK(u != v);
    CHECK(u < n_ && v < n_);
    int ans = numeric_limits<int>::max();
    while (u != v) {
      if (u == -1 || v == -1) return 0; // disconnect
      if (depth_[u] > depth_[v]) swap(u, v);
      ans = min(ans, parent_weight_[v].second);
      v = parent_weight_[v].first;
    }
    return ans;
  }
private:
  int n_;
  vector<pair<V, int>> parent_weight_;
  vector<int> depth_;
};

class tree_minimum_querys {
  using uvcost_t = tuple<V, V, int>;

  vector<uvcost_t> remove_unfocused_edges(vector<uvcost_t>&& random_contraction_edges) {
    //辺を重い順にsort
    sort(random_contraction_edges.begin(), random_contraction_edges.end(), [](const uvcost_t& l, const uvcost_t& r) {
      return get<2>(l) > get<2>(r);
    });

    //最大indexをunion findで使う頂点数に使う
    auto get_max_vertex_index = [](vector<uvcost_t>& v) {
      int ret = 0;
      for (auto& e : v) ret = max(ret, max(get<0>(e), get<1>(e)));
      return ret;
    };
    const int num_vertices = get_max_vertex_index(random_contraction_edges) + 1;

    struct weight_t {
      int w;
      weight_t() : w(0) {}
      weight_t(int w) : w(w) {}
      bool operator<(const weight_t& r) const {
        return w < r.w;
      }
      weight_t& operator+=(const weight_t& r) {
        w = max(w, r.w);
        return *this;
      }
      weight_t operator+(const weight_t& r) const {
        weight_t ret = *this;
        return ret += r;
      }
    };

    weighted_union_find<weight_t> wuf(num_vertices);
    //頂点を重くして、可能な限り焦点を当てている頂点がrootになるように
    FOR(i, focused_num_) wuf.add_weight(i, weight_t(1));

    vector<uvcost_t> ret;
    for (auto& e : random_contraction_edges) {
      V u, v; int w; tie(u, v, w) = e;
      u = wuf.root(u); v = wuf.root(v);
      // if (u == v) continue; 木なので、既に繋がっている辺はこない
      wuf.unite(u, v);
      if (u < focused_num_ && v < focused_num_) {
        // 新しいgomory-hu treeの辺
        // 辺は重い順にsortされているので、この辺で新しく頂点対が繋がった
        //  = gomory_hu treeのedgeの重さは w
        ret.emplace_back(u, v, w);
      }
    }

    return ret;
  }

public:
  // 注目する頂点数
  // random contracionで新しく追加した頂点に対するクエリは飛んでこないので無視出来る
  tree_minimum_querys(int focused_num) : focused_num_(focused_num) {}

  void add_tree(vector<tuple<V, V, int>> random_contraction_edges) {
    //仮の頂点を削除
    random_contraction_edges = remove_unfocused_edges(std::move(random_contraction_edges));
    tmqs_.emplace_back(random_contraction_edges, focused_num_);
  }

  int query(int u, int v) const {
    int ans = numeric_limits<int>::max();
    for (const auto& tmq : tmqs_) {
      ans = min(ans, tmq.query(u, v));
    }
    return ans;
  }
private:
  const int focused_num_;
  vector<tree_minimum_query> tmqs_;
};

class random_contraction {
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

  void contraction(vector<int>& ancestor, vector<unordered_map<V, int>>& contraction_graph_edges, V u, V v, queue<int>& degree_eq_1) {
    assert(!uf_.is_same(u, v));
    const V new_vertex = num_vertices_ + sz(binary_tree_edges_) / 2;
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
    binary_tree_edges_.emplace_back(ua, new_vertex, uw.weight);
    binary_tree_edges_.emplace_back(va, new_vertex, vw.weight);

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
    uf_(g.num_vertices()) {

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
    assert(sz(binary_tree_edges_) <= 2 * (num_vertices_ - 1));
  }

  tree_minimum_querys build_solver() const {
    tree_minimum_querys solver(num_vertices_);
    solver.add_tree(binary_tree_edges_);
    return solver;
  }

  //破壊的操作
  vector<tuple<V, V, int>>&& edges() {
    return move(binary_tree_edges_);
  }

private:
  const int num_vertices_;
  weighted_union_find<degree_weight> uf_;
  vector<tuple<V, V, int>> binary_tree_edges_;
};

class min_cut_query {
public:
  min_cut_query(G& g) : tmq_(g.num_vertices()) {
    fprintf(stderr, "initialize mincut-query tree ...\n");
    for (int i = 0; i < FLAGS_solver_iter; i++) {
      if (i % 10 == 0) fprintf(stderr, "trees .. %d/%d\n", i, FLAGS_solver_iter);
      random_contraction rc(g);
      tmq_.add_tree(rc.edges());
    }
    fprintf(stderr, "completed.\n");
  }

  int query(int u, int v) const {
    int ans = tmq_.query(u, v);
    return ans;
  }

private:
  tree_minimum_querys tmq_;
};
