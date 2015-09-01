#include "generator.h"
#include <set>
#include "base/base.h"
using namespace std;

namespace agl {
unweighted_edge_list gen_erdos_renyi(V num_vertices, size_t avg_deg) {
  avg_deg = min(avg_deg, static_cast<size_t>(max(0, num_vertices - 1)));
  set<pair<V, V>> es;
  std::uniform_int_distribution<V> rng(0, num_vertices - 1);
  while (es.size() < num_vertices * avg_deg) {
    V u = rng(agl::random), v = rng(agl::random);
    if (u == v) continue;
    es.insert(make_pair(u, v));
  }
  return vector<pair<V, V>>(es.begin(), es.end());
};

unweighted_edge_list gen_grid(size_t num_rows, size_t num_cols) {
  auto vid = [=](int i, int j) { return i * num_cols + j; };

  vector<pair<V, V>> es;
  for (auto i : make_irange(num_rows)) {
    for (auto j : make_irange(num_cols)) {
      if (j + 1 < num_cols) es.emplace_back(vid(i, j), vid(i, j + 1));
      if (i + 1 < num_rows) es.emplace_back(vid(i, j), vid(i + 1, j));
    }
  }
  return es;
}

unweighted_edge_list force_undirected(const unweighted_edge_list& es) {
  unweighted_edge_list out(es.size() * 2);
  for (auto i : make_irange(es.size())) {
    V u = es[i].first, v = es[i].second;
    out[i * 2 + 0] = make_pair(u, v);
    out[i * 2 + 1] = make_pair(v, u);
  }
  sort(out.begin(), out.end());
  out.erase(unique(out.begin(), out.end()), out.end());
  return out;
}

unweighted_edge_list agl::gen_barbell(V size_clique) {
  unweighted_edge_list out;
  for (V i : make_irange(2)) {
    for (V v : make_irange(size_clique)) {
      for (V u : make_irange(v)) {
        out.emplace_back(i * size_clique + u, i * size_clique + v);
      }
    }
  }
  out.emplace_back(0, size_clique);
  return out;
}
}  // namespace agl
