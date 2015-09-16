#include "edge_centrality.h"
#include <algorithm>
#include <gflags/gflags.h>
#include "shortest_path/dijkstra.h"
#include "visualize/graphviz.h"
using namespace std;

DEFINE_string(edge_betweenness_centrality_algorithm, "sample", "");
DEFINE_int32(edge_betweenness_centrality_num_samples, 1000, "");

namespace agl {
namespace {
edge_centrality_map init_edge_centrality_map(const G &g) {
  edge_centrality_map res;
  for (auto u : g.vertices()) {
    for (const auto &e : g.edges(u)) {
      res[make_pair(u, to(e))] = 0;
    }
  }
  return res;
}
}  // namespace


edge_centrality_map merge_edge_centrality_map_entries_for_undirected_graph
(const edge_centrality_map &ecm) {
  edge_centrality_map t;

  for (const auto &i : ecm) {
    V u = i.first.first, v = i.first.second;
    if (t.count({v, u})) {
      double &x = t[{v, u}];
      x = (x + i.second) / 2.0;
    } else {
      t[{u, v}] = i.second;
    }
  }
  return t;
}

edge_centrality_map edge_betweenness_centrality(const G &g) {
  if (FLAGS_edge_betweenness_centrality_algorithm == "naive") {
    return edge_betweenness_centrality_naive(g);
  } else if (FLAGS_edge_betweenness_centrality_algorithm == "sample") {
    return edge_betweenness_centrality_sample(g);
  } else if (FLAGS_edge_betweenness_centrality_algorithm == "sample_slow") {
    return edge_betweenness_centrality_sample_slow(g);
  } else {
    CHECK(!"unknown algorithm specified for edge betweenness centrality");
    return init_edge_centrality_map(g);
  }
}

edge_centrality_map edge_betweenness_centrality_naive(const G &g) {
  vector<vector<pair<W, double>>> mat[2];
  for (auto d : directions()) {
    mat[d].resize(g.num_vertices());
    for (auto s : g.vertices()) {
      mat[d][s] = single_source_distance_with_num_paths(g, s, D(d));
    }
  }

  auto res = init_edge_centrality_map(g);
  for (auto u : g.vertices()) {
    for (const auto &e : g.edges(u)) {
      auto v = to(e);
      double sum = 0.0;

      for (auto s : g.vertices()) {
        for (auto t : g.vertices()) {
          const auto &ms = mat[kFwd][s], &mt = mat[kBwd][t];
          assert(is_eq(ms[t].first, mt[s].first));
          assert(is_eq(ms[t].second, mt[s].second));

          const auto &psu = ms[u], &pvt = mt[v], &pst = ms[t];
          if (pst.first == kInfW || psu.first == kInfW || pvt.first == kInfW) continue;
          const W td = psu.first + weight(e) + pvt.first;
          assert(is_le(pst.first, td));
          if (is_lt(pst.first, td)) continue;

          const double np = psu.second * pvt.second / pst.second;
          assert(is_le(np, 1.0));
          sum += np;
        }
      }

      const double ebc = sum / (g.num_vertices() * g.num_vertices());
      res[make_pair(u, v)] = ebc;
    }
  }

  return res;
}


edge_centrality_map edge_betweenness_centrality_sample_slow(const G &g) {
  const auto num_samples = FLAGS_edge_betweenness_centrality_num_samples;
  random_type rng;
  auto res = init_edge_centrality_map(g);

  for (int iter = 0; iter < num_samples; ++iter) {
    const V s = rng(g.num_vertices());
    const V t = rng(g.num_vertices());

    auto ms = single_source_distance_with_num_paths(g, s, kFwd);
    auto mt = single_source_distance_with_num_paths(g, t, kBwd);
    assert(is_eq(ms[t].first, mt[s].first));
    assert(is_eq(ms[t].second, mt[s].second));

    for (auto u : g.vertices()) {
      for (const auto &e : g.edges(u)) {
        auto v = to(e);

        const auto &psu = ms[u], &pvt = mt[v], &pst = ms[t];
        if (pst.first == kInfW || psu.first == kInfW || pvt.first == kInfW) continue;
        const W td = psu.first + weight(e) + pvt.first;
        assert(is_le(pst.first, td));
        if (is_lt(pst.first, td)) continue;

        const double np = psu.second * pvt.second / pst.second;
        assert(is_le(np, 1.0));

        res[make_pair(u, v)] += np / num_samples;
      }
    }
  }

  return res;
}

edge_centrality_map edge_betweenness_centrality_sample(const G &g) {
  const auto num_samples = FLAGS_edge_betweenness_centrality_num_samples;
  auto res = init_edge_centrality_map(g);
  dijkstra_heap<G> hs[2] = {make_dijkstra_heap(g), make_dijkstra_heap(g)};
  vector<pair<W, double>> ds[2];
  ds[kFwd].assign(g.num_vertices(), {kInfW, 0});
  ds[kBwd].assign(g.num_vertices(), {kInfW, 0});
  random_type rng;

  for (int iter = 0; iter < num_samples; ++iter) {
    // Sampled vertex pair
    const V s = rng(g.num_vertices());
    const V t = rng(g.num_vertices());

    hs[kFwd].decrease(s, 0);
    hs[kBwd].decrease(t, 0);
    ds[kFwd][s] = ds[kBwd][t] = make_pair(0, 1);

    // Bidirectional dijkstra search
    vector<pair<V, E>> relaxed_es[2];
    W d = kInfW;
    for (;;) {
      if (hs[kFwd].empty() || hs[kBwd].empty()) break;
      const W lb = hs[kFwd].top_weight() + hs[kBwd].top_weight();
      if (is_le(d, lb)) break;

      D dir = (hs[kFwd].top_weight() < hs[kBwd].top_weight() ? kFwd : kBwd);
      auto &h = hs[dir];
      const V v = h.top_vertex();
      const W w = h.top_weight();
      h.pop();

      for (const auto &e : g.edges(v, dir)) {
        const V tv = to(e);
        const W tw = w + weight(e);
        h.decrease(tv, tw);
        if (is_lt(tw, ds[dir][tv].first)) ds[dir][tv] = make_pair(tw, 0);
        if (is_le(tw, ds[dir][tv].first)) {
          ds[dir][tv].second += ds[dir][v].second;
          relaxed_es[dir].emplace_back(v, e);
        }

        if (ds[1 - dir][tv].first != kInfW) d = min(d, tw + ds[1 - dir][tv].first);
      }
    }
    if (d == kInfW) goto cleanup;

    // Further propagation on relaxed edges
    for (auto dir : directions()) {
      auto &es = relaxed_es[dir];
      reverse(es.begin(), es.end());
      for (const auto &e : es) {
        const V u = e.first, v = to(e.second);
        auto &dsu = ds[1 - dir][u], &dsv = ds[1 - dir][v];
        if (dsv.first == kInfW) continue;
        const W tw = weight(e.second) + dsv.first;
        if (is_lt(tw, dsu.first)) dsu = make_pair(tw, 0);
        if (is_le(tw, dsu.first)) dsu.second += dsv.second;
      }
    }

    // Aggregation
    for (auto dir : directions()) {
      for (const auto &e : relaxed_es[dir]) {
        V u = e.first, v = to(e.second);
        if (dir == kBwd) swap(u, v);

        const auto &psu = ds[kFwd][u], &pvt = ds[kBwd][v], &pst = ds[kFwd][t];
        if (pst.first == kInfW || psu.first == kInfW || pvt.first == kInfW) continue;
        const W td = psu.first + weight(e.second) + pvt.first;
        assert(is_le(pst.first, td));
        if (is_lt(pst.first, td)) continue;

        const double np = psu.second * pvt.second / pst.second;
        assert(is_le(np, 1.0));
        res[make_pair(u, v)] += np / num_samples;
      }
    }

  cleanup:
    // Reset
    for (auto dir : directions()) {
      for (const auto &e : relaxed_es[dir]) {
        const V u = e.first, v = to(e.second);
        ds[kFwd][u] = ds[kBwd][u] = make_pair(kInfW, 0);
        ds[kFwd][v] = ds[kBwd][v] = make_pair(kInfW, 0);
      }
    }
    ds[kFwd][s] = ds[kBwd][s] = make_pair(kInfW, 0);
    ds[kFwd][t] = ds[kBwd][t] = make_pair(kInfW, 0);

    hs[kFwd].clear();
    hs[kBwd].clear();
  }

  return res;
}
}  // namespace agl
