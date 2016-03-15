#include "easy_cui.h"
#include "picojson.h"
using namespace agl;
using namespace picojson;
DEFINE_string(jlog_file, "", "JLOG file path");

string load_json(string& f) {
  ifstream ifs(f);
  if (ifs.fail()) {
    cerr << "faled" << endl;
    return "";
  }
  int begin = static_cast<int>(ifs.tellg());
  ifs.seekg(0, ifs.end);
  int end = static_cast<int>(ifs.tellg());
  int size = end - begin;
  ifs.clear();
  ifs.seekg(0, ifs.beg);
  char* str = new char[size + 1];
  str[size] = '\0';
  ifs.read(str, size);
  return str;
}

unweighted_edge_list shrink(const G& g, vector<V> centers, W rad) {
  int N = g.num_vertices();
  vector<W> mapper(N), dist(N, N);
  V cnt = 0;
  queue<V> que;
  for (const auto& c : centers) {
    que.push(c);
    dist[c] = 0;
    mapper[c] = cnt;
    cnt++;
  }

  for (W d = 0; d < rad; ++d) {
    int s = que.size();
    while (s--) {
      V v = que.front();
      que.pop();
      if (dist[v] == rad)
        for (const auto& u : g.neighbors(v)) {
          if (dist[u] <= d + 1) continue;
          dist[u] = d + 1;
          mapper[u] = mapper[v];
          que.push(u);
        }
    }
  }

  set<pair<V, V>> s;
  for (const auto e : g.edge_list()) {
    V v = mapper[e.first], u = mapper[e.second];
    if (v == u) continue;
    if (v > u) swap(v, u);
    s.insert({v, u});
  }
  unweighted_edge_list es;
  for (const auto& p : s) es.emplace_back(p.first, p.second);
  return es;
}

int main(int argc, char** argv) {
  G g = easy_cui_init(argc, argv);
  CHECK_MSG(FLAGS_force_undirected, "undirected only!!!");
  JLOG_ADD_OPEN("distribution") {
    JLOG_ADD_OPEN(to_string(0).data()) {
      int N = g.num_vertices();
      vector<V> distribution(N, 0);
      for (int i = 0; i < N; ++i) distribution[g.degree(i)]++;
      while (distribution.back() == 0) distribution.pop_back();
      for (int i = 0; i < distribution.size(); ++i)
        JLOG_PUT(to_string(i).data(), distribution[i]);
    }

    string json = load_json(FLAGS_jlog_file);
    value v;
    string err;
    parse(v, json.begin(), json.end(), &err);
    picojson::array& centers_array =
        v.get<object>()["centers"].get<picojson::array>();
    for (auto it = centers_array.begin(); it != centers_array.end(); it++) {
      object& co = it->get<object>();
      for (auto co_it = co.begin(); co_it != co.end(); co_it++) {
        auto rad = co_it->first;
        JLOG_ADD_OPEN(rad.data()) {
          picojson::array& centers = co[rad].get<picojson::array>();
          vector<V> v;
          for (const auto& c : centers) v.push_back(stoi(c.to_str()));
          auto shrinked_es = shrink(g, v, stoi(rad));
          G shg(shrinked_es);
          int N = shg.num_vertices();
          vector<V> distribution(N, 0);
          for (int i = 0; i < N; ++i) distribution[shg.degree(i)]++;
          while (!distribution.empty() && distribution.back() == 0)
            distribution.pop_back();
          for (int i = 0; i < distribution.size(); ++i)
            JLOG_PUT(to_string(i).data(), distribution[i]);
        }
      }
    }
  }

  return 0;
}
