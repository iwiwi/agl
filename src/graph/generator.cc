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

unweighted_edge_list generate_cycle(V num_vertices) {
  unweighted_edge_list es;
  if (num_vertices < 2) {
    return es;
  }
  for (V v = 0; v + 1 < num_vertices; ++v) {
    es.emplace_back(v, v + 1);
  }
  es.emplace_back(num_vertices - 1, 0);
  return es;
}

/**
 * Generate a random scale-free network by the Barabasi-Albert (BA) model.
 * The degree distribution resulting from the BA model is scale free
 * with power-law coefficient &gamma; = 3.
 * \param initial_num is a number of nodes of the initial connected network.
 * \param final_num is a number of finally generated network.
 */
unweighted_edge_list generate_ba(V final_num, V initial_num) {
  assert(initial_num > 2);
  unweighted_edge_list es;
  for (int v = 0; v < initial_num; ++v) {
    for (int u = 0; u < v; ++u) {
      es.emplace_back(u, v);
    }
  }

  for (int v = initial_num; v < final_num; ++v) {
    set<V> next;
    std::uniform_int_distribution<size_t> rng(0, es.size() - 1);
    while (next.size() < (size_t)initial_num) {
      size_t e = rng(agl::random);
      V u = rng(agl::random) % 2 ? es[e].first : es[e].second;
      next.insert(u);
    }
    for (auto u : next) {
      es.emplace_back(u, v);
    }
  }
  return es;
}

/**
 * Generate a random scale-free network by the Dorogovtsev-Mendes-Samukhin (DMS) model.
 * Depending on K0, the power-law coefficient &gamma; takes values from 2 to &infin;.
 * \param initial_num is a number of nodes of the initial connected network.
 * \param final_num is a number of finally generated network.
 * \param K0 is a constant value. The smaller K0 is, the greater &gamma;.
 */
unweighted_edge_list generate_dms(V final_num, V initial_num, V K0) {
  assert(initial_num > K0);
  unweighted_edge_list es;
  vector<V> vs;
  for (int v = 0; v <= initial_num; ++v) {
    for (int u = 0; u < v; ++u) {
      es.emplace_back(u, v);
    }
    for (int i = 0; i < initial_num - K0; ++i) {
      vs.emplace_back(v);
    }
  }

  for (int v = initial_num + 1; v < final_num; ++v) {
    set<V> next;
    std::uniform_int_distribution<size_t> rng(0, vs.size() - 1);
    while (next.size() < (size_t)initial_num) {
      V u = vs[rng(agl::random)];
      next.insert(u);
    }
    for (auto u : next) {
      es.emplace_back(u, v);
      vs.emplace_back(u);
    }
    for (int i = 0; i < initial_num - K0; ++i) {
      vs.emplace_back(v);
    }
  }

  return es;
}

/**
 * Generate a random scale-free network by the Holme-Kim (HK) model.
 * The degree distribution resulting from the HK model is scale free
 * with power-law coefficient &gamma; = 3.
 * \param initial_num is a number of nodes of the initial connected network.
 * \param final_num is a number of finally generated network.
 * \param P is the probability to perform a triangle formation step. 
 */
unweighted_edge_list generate_hk(V final_num, V initial_num, double P) {
  assert(initial_num > 2 && final_num > initial_num);
  unweighted_edge_list es;
  vector<vector<V>> adj(final_num);
  for (int v = 0; v <= initial_num; ++v) {
    for (int u = 0; u < v; ++u) {
      es.emplace_back(u, v);
      adj[u].emplace_back(v);
      adj[v].emplace_back(u);
    }
  }

  std::uniform_real_distribution<> p_rng(0.0, 1.0);
  for (int v = initial_num + 1; v < final_num; ++v) {
    set<V> next;
    V u = -1;
    while (next.size() < (size_t)initial_num) {
      if (next.size() > 0 && p_rng(agl::random) < P) {
        std::uniform_int_distribution<size_t> adj_rng(0, adj[u].size() - 1);
        u = adj[u][adj_rng(agl::random)];
      } else {
        std::uniform_int_distribution<size_t> rng(0, es.size() - 1);
        size_t e = rng(agl::random);
        u = rng(agl::random) % 2 ? es[e].first : es[e].second;
      }
      next.insert(u);
    }
    for (auto u : next) {
      es.emplace_back(u, v);
      adj[u].emplace_back(v);
      adj[v].emplace_back(u);
    }
  }

  return es;
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
