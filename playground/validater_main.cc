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

#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())

DEFINE_string(method, "generate", "generate, test");
DEFINE_string(validation_data_path, "validate.data", "");

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
    for (int cur_level = 0; !path_found && sz(qs) != 0 && sz(qt) != 0; cur_level++) {
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

  int dfs(int v, int t, bool use_slevel, int f) {
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

  vector<pair<int, int>> level;
  vector<int> iter;
  vector<vector<E>> e;
};

G to_directed_graph(G g) {
  vector<pair<V, V>> ret;
  for (auto& e : g.edge_list()) {
    if (e.first < to(e.second)) ret.emplace_back(e.first, to(e.second));
  }
  return G(ret);
}

void generate(const G& g) {
  FILE* fp = fopen(FLAGS_validation_data_path.c_str(), "w");
  dinic dc_base(g.num_vertices());
  FOR(v, g.num_vertices()) for (auto& e : g.edges(v)) {
    dc_base.add_undirected_edge(v, to(e), 1);
  }

  int counter = 0;
  int all = g.num_vertices() * (g.num_vertices() - 1) / 2;
  FOR(s, g.num_vertices()) for (int t = s + 1; t < g.num_vertices(); t++) {
    auto dc = dc_base;
    if (counter % 100 == 0) {
      fprintf(stderr, "count/all = %d/%d\n", counter, all);
    }
    int flow = dc.max_flow(s, t);
    fprintf(fp, "%d %d %d\n", s, t, flow);
    counter++;
  }
  fclose(fp);
}

int main(int argc, char** argv) {

  G g = easy_cui_init(argc, argv);
  g = to_directed_graph(g);

  if (FLAGS_method == "generate") {
    generate(g);
  } else {
    fprintf(stderr, "unrecognized option '%s'\n", FLAGS_method.c_str());
    exit(-1);
  }
}
