#include <easy_cui.h>

#ifdef _WIN32
#pragma comment(linker, "/STACK:3200000") 
#endif

#ifdef _WIN32
template<typename GraphType = G>
GraphType easy_cui_init(int argc, char **argv) {
  JLOG_INIT(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_type == "auto")
    FLAGS_type = guess_type();

  //agl形式は隣接リストを介すると効率が悪いため、直接GraphTypeを生成する
  if (FLAGS_type == "agl") {
    auto g = read_graph_binary<G>(FLAGS_graph.c_str());
    if (FLAGS_force_undirected) {
      g = GraphType(make_undirected(g.edge_list()));
    }
    pretty_print(g);
    return g;
  }

  G::edge_list_type es;
  if (FLAGS_type == "tsv") {
    es = read_edge_list_tsv(FLAGS_graph.c_str());
  } else if (FLAGS_type == "built_in") {
    es = built_in_edge_list(FLAGS_graph.c_str());
  } else if (FLAGS_type == "gen") {
    istringstream iss(FLAGS_graph);
    string family;
    iss >> family;
    if (family == "barbell") {
      V n;
      if (!(iss >> n)) n = 4;
      es = generate_barbell(n);
    } else if (family == "grid") {
      size_t r, c;
      if (!(iss >> r)) r = 4;
      if (!(iss >> c)) c = r;
      es = generate_grid(r, c);
    } else if (family == "erdos_renyi") {
      V n;
      double d;
      if (!(iss >> n)) n = 10;
      if (!(iss >> d)) d = 3.0;
      es = generate_erdos_renyi(n, d);
    } else if (family == "random_planar") {
      V n;
      size_t e;
      if (!(iss >> n)) n = 10;
      if (!(iss >> e)) e = 25;
      es = generate_random_planar(n, e);
    } else if (family == "cycle") {
      V n;
      if (!(iss >> n)) n = 10;
      es = generate_cycle(n);
    } else if (family == "ba") {
      V n, m;
      if (!(iss >> n)) n = 10;
      if (!(iss >> m)) m = 5;
      es = generate_ba(n, m);
    } else if (family == "dms") {
      V n, m, k0;
      if (!(iss >> n)) n = 10;
      if (!(iss >> m)) m = 5;
      if (!(iss >> k0)) k0 = -2;
      es = generate_dms(n, m, k0);
    } else if (family == "hk") {
      V n, m;
      double p;
      if (!(iss >> n)) n = 10;
      if (!(iss >> m)) m = 5;
      if (!(iss >> p)) p = 0.5;
      es = generate_hk(n, m, p);
    } else if (family == "ws") {
      V n, d;
      double p;
      if (!(iss >> n)) n = 10;
      if (!(iss >> d)) d = 4;
      if (!(iss >> p)) p = 0.5;
      es = generate_ws(n, d, p);
    } else if (family == "kronecker") {
      int scale, n;
      size_t avg_deg;
      if (!(iss >> scale)) scale = 5;
      if (!(iss >> avg_deg)) avg_deg = 16;
      vector<vector<double>> mat;
      if (!(iss >> n)) {
        n = 2;
        mat = vector<vector<double>>(n, vector<double>(n));
        mat[0][0] = 0.57;
        mat[0][1] = 0.19;
        mat[1][0] = 0.19;
        mat[1][1] = 0.05;
      }
      for (int i = 0; i < n * n; ++i) {
        if (!(iss >> mat[i / n][i % n])) mat[i / n][i % n] = 1.0 / (n * n);
      }
      es = generate_kronecker(scale, avg_deg, mat);
    } else if (family == "flower") {
      V required, u, v;
      if (!(iss >> required)) required = 44;
      if (!(iss >> u)) u = 2;
      if (!(iss >> v)) v = 2;
      es = generate_uv_flower(required, u, v);
    } else if (family == "shm") {
      V required_num, initial_num;
      int t;
      if (!(iss >> required_num)) required_num = 101;
      if (!(iss >> initial_num)) initial_num = 5;
      if (!(iss >> t)) t = 2;
      es = generate_shm(required_num, initial_num, t);
    } else {
      FAIL_MSG("Unknown generator family: " + family);
    }
  }

  if (FLAGS_force_undirected) es = make_undirected(es);

  GraphType g(es);
  pretty_print(g);
  return g;
}
#endif

DEFINE_int32(num_query, 1000, "");
DEFINE_int64(node_pair_random_seed, 922337203685477583LL, "");
DEFINE_string(method, "gusfield", "both, gomory_hu, gusfield");

#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())

class dinic {
  struct E {
    int to, rev, cap;
    E(int to, int rev, int cap) : to(to), rev(rev), cap(cap) {}
  };

  bool two_sided_bfs(int s, int t) {
    queue<int> qs, qt;
    qs.push(s); qt.push(t);
    level[s].first = level[t].second = 0;
    bfs_revision[s] = s_side_bfs_revision;
    bfs_revision[t] = t_side_bfs_revision;

    int slevel = 0, tlevel = 0;
    while (sz(qs) != 0 && sz(qt) != 0) {
      if (sz(qs) < sz(qt)) {
        int size = sz(qs);
        FOR(_, size) {
          const int v = qs.front(); qs.pop();
          for (auto t : e[v]) {
            if (t.cap == 0 || bfs_revision[t.to] == s_side_bfs_revision) continue;
            if (bfs_revision[t.to] == t_side_bfs_revision) {
              return true; // path was found.
            }
            bfs_revision[t.to] = s_side_bfs_revision;
            level[t.to].first = slevel + 1;
            qs.push(t.to);
          }
        }
        slevel++;
      } else {
        int size = sz(qt);
        FOR(_, size) {
          const int v = qt.front(); qt.pop();
          for (auto t : e[v]) {
            if (e[t.to][t.rev].cap == 0 || bfs_revision[t.to] == t_side_bfs_revision) continue;
            if (bfs_revision[t.to] == s_side_bfs_revision) {
              return true; // path was found.
            }
            bfs_revision[t.to] = t_side_bfs_revision;
            level[t.to].second = tlevel + 1;
            qt.push(t.to);
          }
        }
      }
      tlevel++;
    }
    return false;
  }

  int dfs(int v, int t, bool use_slevel, int f) {
    if (v == t) return f;
    if (dfs_revision[v] != bfs_revision[v]) {
      dfs_revision[v] = bfs_revision[v];
      iter[v] = 0;
    }
    for (int &i = iter[v]; i < sz(e[v]); i++) {
      E& _e = e[v][i];
      if (_e.cap == 0 || bfs_revision[_e.to] / 2 != s_side_bfs_revision / 2) continue;

      bool rec;
      if (use_slevel) rec = bfs_revision[_e.to] == t_side_bfs_revision || level[v].first < level[_e.to].first;
      else rec = level[v].second > level[_e.to].second;
      if (!rec) continue;

      bool next_slevel = use_slevel && bfs_revision[_e.to] == s_side_bfs_revision;
      int d = dfs(_e.to, t, next_slevel, min(f, _e.cap));
      if (d > 0) {
        _e.cap -= d;
        e[_e.to][_e.rev].cap += d;
        return d;
      }
    }
    return 0;
  }

public:
  dinic(const int n)
    : level(n), iter(n), bfs_revision(n), dfs_revision(n), e(n) {
  }

  void add_undirected_edge(int f, int t, int c) {
    e[f].push_back(E(t, sz(e[t]), c));
    e[t].push_back(E(f, sz(e[f]) - 1, c));
  }

  int max_flow(int s, int t) {
    assert(s != t);
    int flow = 0;
    s_side_bfs_revision = 2; t_side_bfs_revision = 3;
    for (; ; s_side_bfs_revision += 2, t_side_bfs_revision += 2) {
      bool path_found = two_sided_bfs(s, t);
      if (!path_found) return flow;
      int f;
      while ((f = dfs(s, t, true, numeric_limits<int>::max())) > 0) {
        flow += f;
      }
    }
  }

  vector<pair<int, int>> level;
  vector<int> iter;
  vector<int> bfs_revision, dfs_revision;
  vector<vector<E>> e;
  int s_side_bfs_revision, t_side_bfs_revision;
};


class dinic_twosided {
  class E {
    const int init_cap_;
    int cap_, revision_;
  public:
    E(int to, int reverse, int cap) :
      init_cap_(cap), cap_(cap), revision_(0), to(to), reverse(reverse) {
    }

    const int to, reverse;
    const int cap(int currenct_revision) {
      if (revision_ != currenct_revision) {
        revision_ = currenct_revision;
        cap_ = init_cap_;
      }
      return cap_;
    }
    void add_cap(int val, int currenct_revision) {
      revision_ = currenct_revision;
      cap_ += val;
    }

    void reset() {
      revision_ = 0;
      cap_ = init_cap_;
    }
  };

  bool two_sided_bfs(int s, int t) {
    queue<int> qs, qt;
    qs.push(s); qt.push(t);
    level[s].first = level[t].second = 0;
    bfs_revision[s] = s_side_bfs_revision;
    bfs_revision[t] = t_side_bfs_revision;

    int slevel = 0, tlevel = 0;
    while (sz(qs) != 0 && sz(qt) != 0) {
      if (sz(qs) < sz(qt)) {
        int size = sz(qs);
        FOR(_, size) {
          const int v = qs.front(); qs.pop();
          for (auto t : e[v]) {
            if (t.cap(graph_revision) == 0 || bfs_revision[t.to] == s_side_bfs_revision) continue;
            if (bfs_revision[t.to] == t_side_bfs_revision) {
              return true; // path was found.
            }
            bfs_revision[t.to] = s_side_bfs_revision;
            level[t.to].first = slevel + 1;
            qs.push(t.to);
          }
        }
        slevel++;
      } else {
        int size = sz(qt);
        FOR(_, size) {
          const int v = qt.front(); qt.pop();
          for (auto t : e[v]) {
            if (e[t.to][t.reverse].cap(graph_revision) == 0 || bfs_revision[t.to] == t_side_bfs_revision) continue;
            if (bfs_revision[t.to] == s_side_bfs_revision) {
              return true; // path was found.
            }
            bfs_revision[t.to] = t_side_bfs_revision;
            level[t.to].second = tlevel + 1;
            qt.push(t.to);
          }
        }
      }
      tlevel++;
    }
    return false;
  }

  int dfs(int v, int t, bool use_slevel, int f) {
    if (v == t) return f;
    if (dfs_revision[v] != bfs_revision[v]) {
      dfs_revision[v] = bfs_revision[v];
      iter[v] = 0;
    }
    for (int &i = iter[v]; i < sz(e[v]); i++) {
      E& _e = e[v][i];
      if (_e.cap(graph_revision) == 0 || bfs_revision[_e.to] / 2 != s_side_bfs_revision / 2) continue;

      bool rec;
      if (use_slevel) rec = bfs_revision[_e.to] == t_side_bfs_revision || level[v].first < level[_e.to].first;
      else rec = level[v].second > level[_e.to].second;
      if (!rec) continue;

      bool next_slevel = use_slevel && bfs_revision[_e.to] == s_side_bfs_revision;
      int d = dfs(_e.to, t, next_slevel, min(f, _e.cap(graph_revision)));
      if (d > 0) {
        _e.add_cap(-d, graph_revision);
        e[_e.to][_e.reverse].add_cap(d, graph_revision);
        return d;
      }
    }
    return 0;
  }

  void add_undirected_edge(int f, int t, int c) {
    e[f].push_back(E(t, sz(e[t]), c));
    e[t].push_back(E(f, sz(e[f]) - 1, c));
  }

  void reset_rivision() {
    FOR(v, n) for (auto& e_ : e[v]) e_.reset();
    memset(bfs_revision.data(), 0, sizeof(bfs_revision[0]) * bfs_revision.size());
    memset(dfs_revision.data(), 0, sizeof(dfs_revision[0]) * dfs_revision.size());
    s_side_bfs_revision = 2;
    t_side_bfs_revision = 3;
    graph_revision = 0;
  }

public:
  dinic_twosided() : n(0) {}
  dinic_twosided(const G& g)
    : n(g.num_vertices()), level(n), iter(n), bfs_revision(n), dfs_revision(n), e(n),
    s_side_bfs_revision(2), t_side_bfs_revision(3), graph_revision(0) {
    FOR(v, n) for (auto& e : g.edges(v)) {
      add_undirected_edge(v, to(e), 1);
    }
  }
  dinic_twosided(const vector<pair<V, V>>& edges, int num_vs)
    : n(num_vs), level(n), iter(n), bfs_revision(n), dfs_revision(n), e(n),
    s_side_bfs_revision(2), t_side_bfs_revision(3), graph_revision(0) {
    for (auto& uv : edges) {
      add_undirected_edge(uv.first, uv.second, 1);
    }

  }

  int max_flow(int s, int t) {
    assert(s != t);
    int flow = 0;
    s_side_bfs_revision += 2;
    t_side_bfs_revision += 2;
    for (; ; s_side_bfs_revision += 2, t_side_bfs_revision += 2) {
      bool path_found = two_sided_bfs(s, t);
      if (!path_found) return flow;
      int f;
      while ((f = dfs(s, t, true, numeric_limits<int>::max())) > 0) {
        flow += f;
      }
    }
    return flow;
  }

  //フローを流す前に実行する
  void reset_graph() {
    graph_revision++;
    if (s_side_bfs_revision >= numeric_limits<decltype(s_side_bfs_revision)>::max() / 2) {
      reset_rivision();
    }
  }

  const int n;
  vector<pair<int, int>> level;
  vector<int> iter;
  vector<int> bfs_revision, dfs_revision;
  vector<vector<E>> e;
  int s_side_bfs_revision, t_side_bfs_revision;
  int graph_revision;

};


class Gusfield {
  int query_dfs(V v, V t, int cost, V par = -1) const {
    if (v == t) return cost;
    for (const auto& to_cost : binary_tree_edges_[v]) {
      V to; int edge_cost; tie(to, edge_cost) = to_cost;
      if (to == par) continue;
      int ncost = query_dfs(to, t, min(cost, edge_cost), v);
      if (ncost != -1) return ncost;
    }

    return -1;
  }

public:

  Gusfield(G g) : num_vertices_(g.num_vertices()), binary_tree_edges_(g.num_vertices()) {
    union_find uf(num_vertices_);
    FOR(v, num_vertices_) for (auto& e : g.edges(v)) {
      uf.unite(v, to(e));
    }
    vector<int> p(num_vertices_);
    FOR(v, num_vertices_) {
      if (v == uf.root(v)) p[v] = -1;
      else p[v] = uf.root(v);
    }

    map<int, dinic_twosided> mp;
    FOR(v, num_vertices_) {
      if (v == uf.root(v)) {
        vector<pair<V, V>> edges;
        FOR(u, num_vertices_) {
          if (!uf.is_same(u, v)) continue;
          for (auto& e : g.edges(u)) {
            edges.emplace_back(u, to(e));
          }
        }
        mp.insert(make_pair(v, dinic_twosided(edges, num_vertices_)));
      }
    }

    FOR(s, num_vertices_) {
      if (p[s] == -1) continue;
      V t = p[s];
      dinic_twosided dc = mp[uf.root(t)];

      int cost = dc.max_flow(s, t);
      //       fprintf(stderr, "(%d,%d) cost = %d\n", s, t, cost);
      binary_tree_edges_[s].emplace_back(t, cost);
      binary_tree_edges_[t].emplace_back(s, cost);

      vector<char> used(num_vertices_);
      queue<int> q;
      q.push(s);
      used[s] = true;
      while (!q.empty()) {
        V v = q.front(); q.pop();
        for (auto& e : dc.e[v]) {
          if (e.cap(dc.graph_revision) == 0 || used[e.to]) continue;
          used[e.to] = true;
          q.push(e.to);
          if (p[e.to] == t) p[e.to] = s;
        }
      }
    }
  }

  int query(V u, V v) const {
    CHECK(u != v);
    CHECK(u < num_vertices_ && v < num_vertices_);
    int ans = query_dfs(u, v, numeric_limits<int>::max());
    if (ans == -1) {
      return 0; // 到達できなかった
    }
    return ans;
  }

private:
  const int num_vertices_;
  vector<vector<pair<V, int>>> binary_tree_edges_;
};

class OptimizedGusfield {
  void build_depth() {
    depth_[0] = 0;

    stack<int> stk;
    FOR(v, num_vertices_) {
      while (depth_[v] == -1) {
        stk.push(v);
        v = parent_weight_[v].first;
      }
      while (!stk.empty()) {
        int u = stk.top(); stk.pop();
        depth_[u] = depth_[v] + 1;
        v = u;
      }
    }
  }

public:

  OptimizedGusfield(vector<pair<V, V>>& edges, int num_vs) :
    num_vertices_(num_vs),
    parent_weight_(num_vs, make_pair(-1, 0)),
    depth_(num_vs, -1) {
    for (V v = 1; v < num_vertices_; v++) {
      parent_weight_[v].first = 0;
    }
    dinic_twosided dc_base(edges, num_vs);

    vector<int> degree(num_vertices_);
    for (auto& e : edges) degree[e.first]++, degree[e.second]++;

    vector<int> used(num_vertices_);
    FOR(s, num_vertices_) {
      if (parent_weight_[s].first == -1) continue;
      V t = parent_weight_[s].first;

      auto dc = dc_base;
      int cost = dc.max_flow(s, t);
      fprintf(stderr, "(%d,%d) cost = %d\n", s, t, cost);
      parent_weight_[s].second = cost;

      if (degree[s] == cost) continue;

      queue<int> q;
      q.push(s);
      const int F = s + 1;
      used[s] = F;
      while (!q.empty()) {
        V v = q.front(); q.pop();
        for (auto& e : dc.e[v]) {
          if (e.cap(dc.graph_revision) == 0 || used[e.to] == F) continue;
          used[e.to] = F;
          q.push(e.to);
          if (parent_weight_[e.to].first == t) parent_weight_[e.to].first = s;
        }
      }
    }

    build_depth();
  }

  int query(V u, V v) {
    CHECK(u != v);
    CHECK(u < num_vertices_ && v < num_vertices_);
    int ans = numeric_limits<int>::max();
    while (u != v) {
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

private:
  const int num_vertices_;
  vector<pair<V, int>> parent_weight_;
  vector<int> depth_;
};

class OptimizedGusfield2 {
  void build_depth() {
    depth_[root_node_] = 0;

    stack<int> stk;
    FOR(v, num_vertices_) {
      while (depth_[v] == -1) {
        stk.push(v);
        v = parent_weight_[v].first;
      }
      while (!stk.empty()) {
        int u = stk.top(); stk.pop();
        depth_[u] = depth_[v] + 1;
        v = u;
      }
    }
  }

public:

  OptimizedGusfield2(vector<pair<V, V>>& edges, int num_vs) :
    num_vertices_(num_vs),
    parent_weight_(num_vs, make_pair(-1, 0)),
    depth_(num_vs, -1) {

    vector<int> degree(num_vertices_);
    for (auto& e : edges) degree[e.first]++, degree[e.second]++;

    root_node_ = -1;
    FOR(v, num_vertices_) if (degree[v] != 2) {
      root_node_ = v; break;
    }
    if (root_node_ == -1) root_node_ = 0;

    vector<int> mincut_order;
    FOR(v, num_vertices_) {
      if (v == root_node_) continue;
      parent_weight_[v].first = root_node_;
      if (degree[v] == 2) {
        parent_weight_[v].second = 2;
      } else {
        mincut_order.push_back(v);
      }
    }
    parent_weight_[root_node_].first = -1;

    dinic_twosided dc_base(edges, num_vs);

    vector<int> used(num_vertices_);
    for (V s : mincut_order) {
      V t = parent_weight_[s].first;

      auto dc = dc_base;
      int cost = dc.max_flow(s, t);
      fprintf(stderr, "(%d,%d) cost = %d\n", s, t, cost);
      parent_weight_[s].second = cost;

      if (degree[s] == cost) continue;

      queue<int> q;
      q.push(s);
      const int F = s + 1;
      used[s] = F;
      while (!q.empty()) {
        V v = q.front(); q.pop();
        for (auto& e : dc.e[v]) {
          if (e.cap(dc.graph_revision) == 0 || used[e.to] == F) continue;
          used[e.to] = F;
          q.push(e.to);
          if (parent_weight_[e.to].first == t) parent_weight_[e.to].first = s;
        }
      }
    }

    build_depth();
  }

  int query(V u, V v) {
    CHECK(u != v);
    CHECK(u < num_vertices_ && v < num_vertices_);
    int ans = numeric_limits<int>::max();
    while (u != v) {
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

private:
  const int num_vertices_;
  vector<pair<V, int>> parent_weight_;
  vector<int> depth_;

  int root_node_;
};

template<class handler_t>
class ConnectedComponentsFilter {
public:
  ConnectedComponentsFilter(const G& g) : n(g.num_vertices()), uf(n), local_indices(n), handlers_indices(n) {

    FOR(v, n) for (auto e : g.edges(v)) {
      V u = to(e);
      uf.unite(u, v);
    }
    int num_connected_components = 0;
    FOR(v, n) {
      if (uf.root(v) != v) local_indices[v] = ++local_indices[uf.root(v)];
      else handlers_indices[v] = num_connected_components++;
    }
    fprintf(stderr, "num_connected_components = %d\n", num_connected_components);

    vector<bool> used(n);
    FOR(v, n) {
      if (uf.root(v) != v) continue;
      used[v] = true;

      const int num_vs = local_indices[uf.root(v)] + 1;
      local_indices[uf.root(v)] = 0;
      vector<pair<V, V>> edges;
      queue<int> q;
      q.push(v);
      while (!q.empty()) {
        V u = q.front(); q.pop();
        FOR(dir, 2) for (auto& e : g.edges(u, D(dir))) {
          V w = to(e);
          if (!used[w]) {
            used[w] = true;
            q.push(w);
          }
          if (dir == 0) {
            edges.emplace_back(local_indices[u], local_indices[w]);
          }
        }
      }

      fprintf(stderr, "root = %d num_vs = %d, edge_size = %d\n", v, num_vs, sz(edges));
      handlers.emplace_back(edges, num_vs);
    }
  }

  int query(V u, V v) {
    if (!uf.is_same(u, v)) return 0;
    int lu = local_indices[u], lv = local_indices[v];
    auto& handler = handlers[handlers_indices[uf.root(u)]];
    if (lu == lv) {
      printf("u = %d, lu = %d, v = %d, lv = %d\n", u, lu, v, lv);
      CHECK(false);
    }
    return handler.query(lu, lv);
  }

private:
  int n;
  union_find uf;
  vector<int> local_indices, handlers_indices;
  vector<handler_t> handlers;
};


// 2ECC = 2-edge connected components
class OptimizedGusfieldWith2ECC {
public:
  const int n;
  const G& g;

  void lowlink_dfs(int v, int par, int& cur_ord) {
    lowlink[v] = order[v] = cur_ord++;
    FOR(dir, 2) for (auto to : g.edges(v, D(dir))) {
      if (to == par) continue;
      if (order[to] == -1) {
        lowlink_dfs(to, v, cur_ord);
        lowlink[v] = min(lowlink[v], lowlink[to]);
        if (order[v] < lowlink[to]) bridge.emplace_back(v, to);
        else biconnected_graphs_edges.emplace_back(v, to);
      } else {
        lowlink[v] = min(lowlink[v], lowlink[to]);
        if (v < to) biconnected_graphs_edges.emplace_back(v, to);
      }
    }
  }

  OptimizedGusfieldWith2ECC(const G& g) : n(g.num_vertices()), g(g), uf(n), lowlink(n, -1), order(n, -1) {
    FOR(v, n) for (auto& e : g.edges(v)) {
      uf.unite(v, to(e));
    }

    const int num_edges = g.num_edges();

    FOR(v, n) if (uf.root(v) == v) {
      int cur_ord = 0;
      lowlink_dfs(v, -1, cur_ord);
    }

    CHECK(sz(bridge) + sz(biconnected_graphs_edges) == num_edges);

    G new_g(biconnected_graphs_edges, n);
    biconnected_graph_handler.reset(new ConnectedComponentsFilter<OptimizedGusfield2>(new_g));
  }

public:
  int query(V u, V v) {
    int ans = biconnected_graph_handler->query(u, v);
    if (ans == 0) {
      if (uf.is_same(u, v)) return 1; // 橋で間接的につながっている
      else return 0;
    }
    return ans;
  }

private:
  union_find uf;
  vector<int> lowlink, order;
  vector<pair<V, V>> bridge, biconnected_graphs_edges;

  unique_ptr<ConnectedComponentsFilter<OptimizedGusfield2>> biconnected_graph_handler;
};

G to_directed_graph(G g) {
  vector<pair<V, V>> ret;
  for (auto& e : g.edge_list()) {
    if (e.first < to(e.second)) ret.emplace_back(e.first, to(e.second));
  }
  return G(ret);
}

using Gusfield2 = ConnectedComponentsFilter<OptimizedGusfield>;
using Gusfield3 = OptimizedGusfieldWith2ECC;

bool check_min_cut_query(Gusfield2& gf2, Gusfield3& gf3, int s, int t, const G& g) {
  // int naive_w = dc.max_flow(s, t);
  int gf2_w = gf2.query(s, t);
  int gf3_w = gf3.query(s, t);

  if (gf2_w != gf3_w) {
    fprintf(stderr, "unmatched. (S,T) = (%d,%d), naive = None, gf2_w = %d, gf3_w = %d\n", s, t, gf2_w, gf3_w);
  }
  return gf2_w == gf3_w;
}

void test(const G& g) {
  Gusfield3 gf3(g);
  Gusfield2 gf2(g);
  xorshift64star gen_node(FLAGS_node_pair_random_seed);

  if (g.num_vertices() > 5000) {
    int unmatch = 0;
    for (int counter = 0; counter < FLAGS_num_query; counter++) {
      V s = gen_node() % g.num_vertices();
      V t = gen_node() % (g.num_vertices() - 1);
      if (s <= t) t++;
      if (counter % 100 == 0) {
        fprintf(stderr, "count/unmatch/all : %d/%d/%d, \n", counter, unmatch, FLAGS_num_query);
      }
      bool is_matched = check_min_cut_query(gf2, gf3, s, t, g);
      if (!is_matched) unmatch++;
    }
    JLOG_PUT("result.all", FLAGS_num_query);
    JLOG_PUT("result.match", (FLAGS_num_query - unmatch));
    JLOG_PUT("result.unmatch", unmatch);
  } else {
    int counter = 0;
    int unmatch = 0;
    FOR(s, g.num_vertices()) for (int t = s + 1; t < g.num_vertices(); t++) {
      if (counter % 100 == 0) {
        fprintf(stderr, "count/unmatch/all : %d/%d/%d, \n", counter, unmatch, -1);
      }
      bool is_matched = check_min_cut_query(gf2, gf3, s, t, g);
      if (!is_matched) unmatch++;
      counter++;
    }

    JLOG_PUT("result.all", counter);
    JLOG_PUT("result.match", (counter - unmatch));
    JLOG_PUT("result.unmatch", unmatch);
  }
}

void tester() {
  test(to_directed_graph(built_in_graph("ca_grqc")));
  exit(0);
}

int main(int argc, char** argv) {
  // tester();

  G g = easy_cui_init(argc, argv);
  g = to_directed_graph(g);

  if (FLAGS_method == "test") {
    test(g);
  } else if (FLAGS_method == "gusfield") {
    Gusfield3 gf(g);
  } else {
    fprintf(stderr, "unrecognized option '%s'\n", FLAGS_method.c_str());
    exit(-1);
  }
}
