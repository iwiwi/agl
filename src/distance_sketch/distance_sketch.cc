#include "distance_sketch.h"
#include "agl.h"
using namespace std;

namespace agl {
namespace distance_sketch {
///////////////////////////////////////////////////////////////////////////////
// Rank
///////////////////////////////////////////////////////////////////////////////
vector<rank_type> generate_rank_array(V num_vertices) {
  vector<rank_type> rank(num_vertices);
  for (V v = 0; v < num_vertices; ++v) rank[v] = agl::random();
  return rank;
}

void pretty_print_rank_array(const std::vector<rank_type>& rank_array, ostream &ofs) {
  ofs << "[";
  for (V v : make_irange(rank_array.size())) {
    auto r = rank_array[v];
    ofs << v << ":" << fixed << setprecision(2) << r / (double)numeric_limits<rank_type>::max() << ", ";
  }
  ofs << "]" << endl;
}


///////////////////////////////////////////////////////////////////////////////
// Static Sketch
///////////////////////////////////////////////////////////////////////////////

vertex_sketch_raw purify_sketch(vertex_sketch_raw sketch, size_t k, const rank_array &ranks) {
  sort(sketch.begin(), sketch.end());

  vertex_sketch_raw new_sketch;
  priority_queue<rank_type> thr;
  for (const auto &e : sketch) {
    rank_type r = ranks[e.v];
    if (thr.size() < k || thr.top() > r) {
      new_sketch.emplace_back(e);
      thr.push(r);
      if (thr.size() > k) thr.pop();
    }
  }

  return new_sketch;
}

vertex_sketch_raw compute_all_distances_sketch_from
(const G& g, V v, size_t k, const std::vector<rank_type> &ranks, D d) {
  vertex_sketch_raw sketch;

  auto f = [&](V v, W w) -> bool {
    sketch.emplace_back(v, w);
    return true;
  };
  visit_by_distance(g, v, f, d);

  return purify_sketch(sketch, k, ranks);
}

all_distances_sketches compute_all_distances_sketches
(const G& g, size_t k, const std::vector<rank_type> &ranks, D d) {
  all_distances_sketches ads;
  ads.k = k;
  ads.ranks = ranks;
  if (ads.ranks.empty()) {
    ads.ranks = generate_rank_array(g.num_vertices());
  }
  assert((V)ads.ranks.size() == g.num_vertices());
  ads.sketches.resize(g.num_vertices());

  vector<priority_queue<pair<W, V>>> thr(g.num_vertices());

  vector<V> ord(g.num_vertices());
  iota(ord.begin(), ord.end(), 0);
  sort(ord.begin(), ord.end(), [&](V u, V v) {
    return ads.ranks[u] < ads.ranks[v];
  });

  for (V v_source : ord) {
    auto f = [&](V v_visiting, W d) -> bool {
      vertex_sketch_raw &sketch = ads.sketches[v_visiting];
      if (!sketch.empty() && sketch.back().v == v_source) return false;
      if (thr[v_visiting].size() >= k && make_pair(d, v_source) > thr[v_visiting].top()) {
        return false;
      }
      sketch.emplace_back(v_source, d);
      thr[v_visiting].push({d, v_source});
      if (thr[v_visiting].size() > k) thr[v_visiting].pop();
      return true;
    };

    visit_by_distance(g, v_source, f, reverse_direction(d));
  }

  for (auto &s : ads.sketches) {
    sort(s.begin(), s.end());
  }
  return ads;
}

void pretty_print(const vertex_sketch_raw &s, std::ostream &ofs) {
  for (const auto e : s) {
    pretty_print(e, ofs);
    ofs << ", ";
  }
  ofs << endl;
}

void pretty_print(const all_distances_sketches& ads, std::ostream& ofs) {
  pretty_print_rank_array(ads.ranks, ofs);
  for (const auto &s : ads.sketches) {
    pretty_print(s, ofs);
  }
}

///////////////////////////////////////////////////////////////////////////////
// Sketch Retrieval Shortcuts
///////////////////////////////////////////////////////////////////////////////
vertex_sketch_raw sketch_retrieval_shortcuts::retrieve_sketch(V v) {
  vertex_sketch_raw ads;
  priority_queue<rank_type> thr;
  auto add_entry = [&](V v_source, W d) -> bool {
      const rank_type r = ranks[v_source];
      if (!ads.empty() && ads.back().v == v_source) return false;
      if (thr.size() >= k && r > thr.top()) return false;
      ads.emplace_back(v_source, d);
      thr.push(r);
      if (thr.size() > k) thr.pop();
      return true;
  };

  // TODO: supporting |std::pair| in |dijkstra_heap|
  typedef pair<W, V> queue_entry_type;
  priority_queue<queue_entry_type, vector<queue_entry_type>, greater<queue_entry_type>> que;
  unordered_map<V, W> pot;  // TODO

  que.emplace(0, v);
  pot[v] = 0;
  while (!que.empty()) {
    V v = que.top().second;
    W d = que.top().first;
    que.pop();
    if (d > pot[v] || !add_entry(v, d)) continue;

    auto sv = retrieve_shortcuts(v);
    for (const auto &e : sv) {
      V tv = e.v;
      W td = d + e.d;
      auto i = pot.insert(make_pair(tv, td));
      if (!i.second && i.first->second <= td) continue;
      que.emplace(td, tv);
      i.first->second = td;
    }
  }
  sort(ads.begin(), ads.end());
  return ads;
}

sketch_retrieval_shortcuts compute_sketch_retrieval_shortcuts
(const G& g, size_t k, const rank_array& ranks, D d) {
  return compute_sketch_retrieval_shortcuts_via_ads_fast(g, k, ranks, d);
}

sketch_retrieval_shortcuts compute_sketch_retrieval_shortcuts_via_ads_naive
(const G& g, size_t k, const rank_array& ranks, D d) {
  all_distances_sketches ads = compute_all_distances_sketches(g, k, ranks, d);

  sketch_retrieval_shortcuts srs(g);
  srs.k = ads.k;
  srs.ranks = ads.ranks;
  srs.sketches.resize(g.num_vertices());

  // (Distance, to, from)
  using tuple_type = tuple<W, V, V>;
  vector<tuple_type> es;
  {
    for (V v : g.vertices()) {
      const vertex_sketch_raw &s = ads.retrieve_sketch(v);
      for (const auto &e : s) {
        es.emplace_back(make_tuple(e.d, e.v, v));
      }
    }
  }
  sort(es.begin(), es.end());

  for (const tuple_type &e : es) {
    const W d = get<0>(e);
    const V s = get<1>(e);
    const V v = get<2>(e);

    vertex_sketch_raw a = srs.retrieve_sketch(v);
    cout << make_tuple(v, s, d) << " " << endl;
    pretty_print(a);
    if (find(a.begin(), a.end(), entry(s, d)) != a.end()) {
      cout << endl;
      continue;
    }
    srs.sketches[v].emplace_back(s, d);
  }

  return srs;
}

sketch_retrieval_shortcuts compute_sketch_retrieval_shortcuts_via_ads_fast
(const G& g, size_t k, const rank_array& ranks, D d) {
  all_distances_sketches ads = compute_all_distances_sketches(g, k, ranks, d);

  sketch_retrieval_shortcuts srs(g);
  srs.k = ads.k;
  srs.ranks = ads.ranks;
  srs.sketches.resize(g.num_vertices());

  // (Distance, to, from)
  using tuple_type = tuple<W, V, V>;
  vector<tuple_type> es;
  vector<vertex_sketch_raw> rev(g.num_vertices());
  {
    for (V v : g.vertices()) {
      const vertex_sketch_raw &s = ads.retrieve_sketch(v);
      for (const auto &e : s) {
        es.emplace_back(make_tuple(e.d, e.v, v));
        rev[e.v].emplace_back(v, e.d);
      }
    }
  }
  sort(es.begin(), es.end());

  priority_queue<tuple_type, vector<tuple_type>, greater<tuple_type>> gs;
  for (const tuple_type &e : es) {
    while (!gs.empty() && gs.top() < e) gs.pop();
    if (!gs.empty() && gs.top() == e) continue;

    const W d = get<0>(e);
    const V s = get<1>(e);
    const V v = get<2>(e);
    if (d == 0) continue;

    srs.sketches[v].emplace_back(s, d);

    for (const auto &r : rev[v]) {
      gs.emplace(make_tuple(r.d + d, s, r.v));  // r.v -> v -> s
    }
  }
  return srs;
}
}  // namespace distance_sketch
}  // namespace agl
