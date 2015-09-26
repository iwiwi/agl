#include "generator.h"
#include <set>
#include "base/base.h"
using namespace std;

namespace agl {
unweighted_edge_list generate_path(V num_vertices) {
  unweighted_edge_list es;
  for (V v = 0; v + 1 < num_vertices; ++v) {
    es.emplace_back(v, v + 1);
  }
  return es;
}

unweighted_edge_list generate_erdos_renyi(V num_vertices, double avg_deg) {
  avg_deg = min(avg_deg, static_cast<double>(max(0, num_vertices - 1)));
  set<pair<V, V>> es;
  std::uniform_int_distribution<V> rng(0, num_vertices - 1);
  while (es.size() < num_vertices * avg_deg) {
    V u = rng(agl::random), v = rng(agl::random);
    if (u == v) continue;
    es.insert(make_pair(u, v));
  }
  return vector<pair<V, V>>(es.begin(), es.end());
};

unweighted_edge_list generate_grid(size_t num_rows, size_t num_cols) {
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

unweighted_edge_list generate_barbell(V size_clique) {
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

unweighted_edge_list generate_random_planar(V num_vertices, size_t num_edges) {
  using namespace agl::geometry2d;

  uniform_real_distribution<double> urd(0.0, 1.0);
  vector<point_type> points(num_vertices);
  for (auto &p : points) p = point_type(urd(random), urd(random));

  vector<tuple<double, V, V>> es_sorted;
  for (V v : make_irange(num_vertices)) {
    for (V u : make_irange(v)) {
      es_sorted.emplace_back(abs(points[u] - points[v]), u, v);
    }
  }
  sort(es_sorted.begin(), es_sorted.end());

  set<pair<V, V>> es;
  for (const auto &t : es_sorted) {
    if (es.size() >= num_edges) break;
    V u = get<1>(t);;
    V v = get<2>(t);
    if (u == v) continue;

    bool ins = false;
    for (const auto &e : es) {
      if (u == e.first || u == e.second ||
          v == e.first || v == e.second) continue;
      ins = does_intersect(segment_type(points[u], points[v]),
                            segment_type(points[e.first], points[e.second]));
      if (ins) break;
    }

    if (!ins) es.emplace(u, v);
    cout << es.size() << endl;
  }

  return unweighted_edge_list(es.begin(), es.end());
}

unweighted_edge_list gen_random_planar_(V num_vertices, size_t num_edges) {
  using namespace agl::geometry2d;

  uniform_real_distribution<double> urd(0.0, 1.0);
  vector<point_type> points(num_vertices);
  for (auto &p : points) p = point_type(urd(random), urd(random));

  set<pair<V, V>> es;
  while (es.size() < num_edges) {
    V u = random(num_vertices);
    V v = random(num_vertices);
    if (u == v) continue;

    bool ins = false;
    for (const auto &e : es) {
      if (u == e.first || u == e.second ||
          v == e.first || v == e.second) continue;
      ins = does_intersect(segment_type(points[u], points[v]),
                            segment_type(points[e.first], points[e.second]));
      if (ins) break;
    }

    if (!ins) es.emplace(u, v);
    cout << es.size() << endl;
  }

  return unweighted_edge_list(es.begin(), es.end());
}

unweighted_edge_list generate_random_spanning_tree(V num_vertices) {
  union_find uf(num_vertices);
  unweighted_edge_list es;
  std::uniform_int_distribution<V> rng(0, num_vertices - 1);

  while (es.size() + 1 < (size_t)num_vertices) {
    V u = rng(agl::random), v = rng(agl::random);
    if (uf.is_same(u, v)) continue;
    es.emplace_back(u, v);
    uf.unite(u, v);
  }

  return es;
}

unweighted_edge_list make_undirected(const unweighted_edge_list& es) {
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

template<>
unweighted_edge_list add_random_weight<unweighted_graph>(const unweighted_edge_list &es) {
  return es;
}

template<>
unweighted_edge_list add_unit_weight<unweighted_graph>(const unweighted_edge_list &es) {
  return es;
}
}  // namespace agl
