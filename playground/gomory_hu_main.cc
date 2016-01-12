#include <easy_cui.h>

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
    for (auto& p : level) p.first = p.second = -1;

    queue<int> qs, qt;
    qs.push(s); qt.push(t);
    level[s].first = level[t].second = 0;

    bool path_found = false;
    for (int cur_level = 0; !path_found && sz(qs) != 0 && sz(qt) != 0 ; cur_level++) {
      {
        int size = sz(qs);
        FOR(_, size) {
          const int v = qs.front(); qs.pop();
          for (auto t : e[v]) {
            if (t.cap == 0) continue;
            int& nlev = level[t.to].first;
            if (nlev != -1) continue;
            int& opplev = level[t.to].second;
            if (opplev == -1) {
              nlev = cur_level + 1;
              qs.push(t.to);
            } else {
              path_found = true;
              break;
            }
          }
        }
      }
      if (path_found) break;
      {
        int size = sz(qt);
        FOR(_, size) {
          const int v = qt.front(); qt.pop();
          for (auto t : e[v]) {
            if (e[t.to][t.rev].cap == 0) continue;
            int& nlev = level[t.to].second;
            if (nlev != -1) continue;
            int& opplev = level[t.to].first;
            if (opplev == -1) {
              nlev = cur_level + 1;
              qt.push(t.to);
            } else {
              path_found = true;
              break;
            }
          }
        }
      }
    }
    return path_found;
  }

  int dfs(int v, int t,bool use_slevel, int f) {
    if (v == t) return f;
    for (int &i = iter[v]; i < sz(e[v]); i++) {
      E& _e = e[v][i];
      if (_e.cap == 0) continue;
      
      bool rec;
      if (use_slevel) rec = level[v].first < level[_e.to].first || level[_e.to].first == -1 && level[_e.to].second != -1;
      else rec = level[v].second > level[_e.to].second;
      if (!rec) continue;

      bool next_slevel = use_slevel && level[v].first < level[_e.to].first;
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
    : level(n), iter(n), e(n) {
  }

  void add_undirected_edge(int f, int t, int c) {
    e[f].push_back(E(t, sz(e[t]), c));
    e[t].push_back(E(f, sz(e[f]) - 1, c));
  }

  int max_flow(int s, int t) {
    assert(s != t);
    int flow = 0;
    while (true) {
      bool path_found = two_sided_bfs(s, t);
      if (!path_found) return flow;
      iter.assign(iter.size(), 0);
      int f;
      while ((f = dfs(s, t, true, numeric_limits<int>::max())) > 0) {
        flow += f;
      }
    }
  }

  vector<pair<int,int>> level;
  vector<int> iter;
  vector<vector<E>> e;
};

class Gomory_Hu {
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

  void rec_flow(vector<tuple<int,int,int>>& g, vector<char> contracted) {

    vector<char> svs_contracted, tvs_contracted;
    vector<tuple<int, int, int>> svs_g, tvs_g;

    {
      vector<char> has(num_vertices_);
      V s = -1, t = -1;
      auto add_sv_pair = [&s, &t](int x) {
        if (s == -1) s = x;
        else if (t == -1 && s != x) t = x;
      };
      for (auto& e : g) {
        V u, v; tie(u, v, std::ignore) = e;
        has[u] = has[v] = true;
        if (!contracted[u]) add_sv_pair(u);
        if (!contracted[v]) add_sv_pair(v);
      }
      if (t == -1) return;

      dinic dc(num_vertices_);
      for (auto& e : g) {
        V u, v; int c; tie(u, v, c) = e;
        dc.add_undirected_edge(u, v, c);
      }

      int cost = dc.max_flow(s, t);
      binary_tree_edges_[s].emplace_back(t, cost);
      binary_tree_edges_[t].emplace_back(s, cost);
      fprintf(stderr, "(%d,%d) cost = %d\n", s, t, cost);

      queue<V> q;
      q.push(s);
      has[s] = false;
      while (!q.empty()) {
        V v = q.front(); q.pop();
        for (auto e : dc.e[v]) {
          if (e.cap == 0 || !has[e.to]) continue;
          q.push(e.to);
          has[e.to] = false;
        }
      }
      
      auto new_graph = [&, n = num_vertices_](const vector<tuple<int,int,int>>& g, bool left) {
        union_find uf(n);
        vector<char> new_contracted = contracted;
        int cut_cost = 0;
        for (auto& e : g) {
          V u, v; int c; tie(u, v, c) = e;
          if (has[u] == has[v] && has[u] == left) {
            uf.unite(u, v);
          }
          if (has[u] != has[v]) {
            cut_cost += c;
          }
          if (has[u] == left) new_contracted[u] = true;
          if (has[v] == left) new_contracted[v] = true;
        }
        assert(cut_cost == cost);
        map<pair<int, int>, int> edge_cost;
        for (auto& e : g) {
          V u, v; int c; tie(u, v, c) = e;
          u = uf.root(u), v = uf.root(v);
          if (u == v) continue;
          if (u > v) swap(u, v);
          edge_cost[pair<int, int>(u, v)] += c;
        }

        vector<tuple<int, int, int>> ret;
        for (auto& kv : edge_cost) {
          ret.emplace_back(kv.first.first, kv.first.second, kv.second);
        }
      
        return make_pair(new_contracted, ret);
      };

      tie(svs_contracted, svs_g) = new_graph(g, true);
      tie(tvs_contracted, tvs_g) = new_graph(g, false);
    }

    // dinicで使うメモリを開放する
    rec_flow(svs_g, svs_contracted);
    rec_flow(tvs_g, tvs_contracted);
  }

public:
  // g の辺を使い、グラフを作成する(勝手にundirectedとして読み替えている)
  Gomory_Hu(G& g) :
    num_vertices_(g.num_vertices()),
    binary_tree_edges_(g.num_vertices()) {
    vector<unordered_set<pair<V, V>>> contraction_graph_edges(g.num_vertices());
    typename G::edge_list_type initial_edges(g.edge_list());
    vector<int> ancestor(g.num_vertices());
    FOR(i, num_vertices_) ancestor[i] = i;

    union_find uf_(num_vertices_);
    for (auto edge : initial_edges) {
      V u = edge.first, v = to(edge.second);
      contraction_graph_edges[u].emplace(u, v);
      contraction_graph_edges[v].emplace(v, u);
      uf_.unite(u, v);
    }
    map<uint32_t, vector<int>> mp;
    FOR(i, num_vertices_) mp[uf_.root(i)].push_back(i);
    for (auto& kv : mp) {
      vector<tuple<int, int, int>> g;
      const int uf_root = uf_.root(kv.second[0]);
      for (auto& e : initial_edges) {
        if (uf_.root(e.first) == uf_root) g.emplace_back(e.first, e.second, 1);
      }
      vector<char> contracted(num_vertices_);
      rec_flow(g, contracted);
    }
    
    
    // gが連結とは限らないので、uv_costs.size() == g.num_vertices() - 1ではない
    assert(int(binary_tree_edges_.size() - g.num_vertices()) <= g.num_vertices() - 1);
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

class Gusfield{
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
  
  Gusfield(G& g) : num_vertices_(g.num_vertices()), binary_tree_edges_(g.num_vertices()) {
    union_find uf(num_vertices_);
    FOR(v, num_vertices_) for (auto& e : g.edges(v)) {
      uf.unite(v, to(e));
    }
    vector<int> p(num_vertices_);
    FOR(v, num_vertices_) {
      if (v == uf.root(v)) p[v] = -1;
      else p[v] = uf.root(v);
    }

    FOR(s, num_vertices_) {
      if (p[s] == -1) continue;
      V t = p[s];
      dinic dc(num_vertices_);
      FOR(v, num_vertices_) {
        if (!uf.is_same(s, v)) continue;
        for (auto& e : g.edges(v)) {
          dc.add_undirected_edge(v, to(e), 1);
        }
      }

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
          if (e.cap == 0 || used[e.to]) continue;
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


G to_directed_graph(G& g) {
  vector<pair<V, V>> ret;
  for (auto& e : g.edge_list()) {
    if (e.first < to(e.second)) ret.emplace_back(e.first, to(e.second));
  }
  return G(ret);
}

void test(G& g) {
  g = to_directed_graph(g);
  puts("--- Gomory_Hu ---");
  Gomory_Hu gh(g);
  puts("--- Gusfield ---");
  Gusfield gf(g);

  dinic dc_base(g.num_vertices());
  for (auto e : g.edge_list()) {
    dc_base.add_undirected_edge(e.first, e.second, 1);
  }

  FOR(v, g.num_vertices()){
    FOR(u, v) {
      auto dc = dc_base;
      int x = dc.max_flow(u, v);
      //int y = gh.query(u, v);
      int y = x;
      int z = gf.query(u, v);
      printf("(%d,%d) dinic : %d, Gomory-Hu : %d, Gusfield : %d\n",u,v,x,y,z);
      assert(x == y && y == z);
    }
  }
}

int main(int argc, char** argv) {
  G g = easy_cui_init(argc, argv);
  g = to_directed_graph(g);

  if(FLAGS_method == "gomory_hu") {
    Gomory_Hu gh(g);
  } else if(FLAGS_method == "gusfield") {
    Gusfield gf(g);
  } else if(FLAGS_method == "both") {
    Gomory_Hu gh(g);
    Gusfield gf(g);
  } else {
    fprintf(stderr, "unrecognized option '%s'\n", FLAGS_method.c_str());
    exit(-1);
  }
}
