#include <easy_cui.h>

using namespace std;

DEFINE_int32(num_query, 10000000, "");
DEFINE_int64(node_pair_random_seed, 922337203685477583LL, "");
DEFINE_string(tree_path, "", "");

class tree_query {

  void build(vector<vector<pair<V, int>>>& edges) {
    depth_.resize(num_vertices_, -1);
    parent_weight_.resize(num_vertices_, make_pair(-2, -2));

    for(V v = 0; v < num_vertices_; v++) {
      if (depth_[v] >= 0) continue;

      depth_[v] = 0;
      parent_weight_[v] = make_pair(-1, -1);
      queue<int> q;
      q.push(v);
      while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto& to_weight : edges[u]) {
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
  static tree_query from_file(const string& path) {
    ifstream ifs(path.c_str());
    CHECK(ifs.is_open());
    return from_file(ifs);
  }

  static tree_query from_file(istream& is) {
    vector<tuple<V, V, int>> input;
    int s, v, weight;
    while (is >> s >> v >> weight) input.emplace_back(s, v, weight);
    return tree_query(std::move(input));
  }
  tree_query() {}

  tree_query(vector<tuple<V, V, int>> input) {
    num_vertices_ = (int)input.size() + 1;

    vector<vector<pair<V, int>>> edges(num_vertices_);
    vector<int> depth_;
    for (auto& e : input) {
      int s, v, weight; tie(s, v, weight) = e;
      edges[s].emplace_back(v, weight);
      edges[v].emplace_back(s, weight);
    }
    input.clear(); input.shrink_to_fit();

    build(edges);
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

  const int num_vertices() { return num_vertices_; }
  
public:
  int num_vertices_;
  vector<pair<V, int>> parent_weight_;
  vector<int> depth_;
};


void from_file() {
  tree_query tq;
  JLOG_PUT_BENCHMARK("initialize_time") {
    tq = tree_query::from_file(FLAGS_tree_path);
  }

  agl::random_type random(FLAGS_node_pair_random_seed);
  JLOG_PUT_BENCHMARK("query_time") {
    for(int i = 0; i < FLAGS_num_query; i++) {
      V s = random() % tq.num_vertices();
      V t = random() % (tq.num_vertices() - 1);
      if(s <= t) t++;
      tq.query(s, t);
    }
  }
}

int main(int argc, char** argv) {
  JLOG_INIT(&argc, argv);
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  from_file();

  return 0;
}
