#include "distance_sketch.h"
#include "agl.h"
using namespace std;

DEFINE_int32(distance_sketch_k, 16, "");
DEFINE_bool(distance_sketch_implicit_neighbor, true, "");

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

W find_distance(const vertex_sketch_raw &sketch, V v ) {
  for (const auto &e : sketch) if (e.v == v) return e.d;
  return kInfW;
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

  visitor_by_distance<G> visitor(g);
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

    visitor.visit(v_source, f, reverse_direction(d));
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
// Estimation
///////////////////////////////////////////////////////////////////////////////

vertex_sketch_raw sort_by_vertices(vertex_sketch_raw sketch) {
  sort(sketch.begin(), sketch.end(),
       [](entry e1, entry e2) { return e1.v < e2.v; });
  return sketch;
}

vector<double> compute_taus
(size_t k, const rank_array &ranks, const vertex_sketch_raw &sketch) {
  const size_t n = sketch.size();
  vector<size_t> ord(n);
  iota(ord.begin(), ord.end(), 0);
  sort(ord.begin(), ord.end(), [&](int i, int j) { return sketch[i] < sketch[j]; });

  vector<double> taus(n);
  priority_queue<double> pq;
  for (size_t k = 0; k < n; ++k) {
    size_t i = ord[k];
    const auto &e = sketch[i];

    taus[i] = pq.size() < k ? 1.0 : pq.top();
    const double p = rank_to_p(ranks[e.v]);
    assert(p <= taus[i]);

    pq.push(p);
    if (pq.size() > k) pq.pop();
  }

  return taus;
}

double estimate_closeness_centrality
(size_t k, const rank_array &ranks, const vertex_sketch_raw &sketch,
 std::function<double(W)> distance_decay_function) {
  double ans = 0.0;
  vector<double> taus = compute_taus(k, ranks, sketch);
  assert(taus.size() == sketch.size());
  for (size_t i = 0; i < sketch.size(); ++i) {
    const auto &e = sketch[i];
    ans += 1.0 / taus[i] * distance_decay_function(e.d);
  }
  return ans;
}

double estimate_distance(const vertex_sketch_raw& sketch1,
                         const vertex_sketch_raw& sketch2) {
  const vertex_sketch_raw s1 = sort_by_vertices(sketch1);
  const vertex_sketch_raw s2 = sort_by_vertices(sketch2);

  size_t i1 = 0, i2 = 0;
  W d = kInfW;
  while (i1 < s1.size() && i2 < s2.size()) {
    const V v1 = s1[i1].v;
    const V v2 = s2[i2].v;
    /**/ if (v1 < v2) ++i1;
    else if (v1 > v2) ++i2;
    else {
      d = min(d, s1[i1].d + s2[i2].d);
      ++i1;
      ++i2;
    }
  }

  return d;
}

///////////////////////////////////////////////////////////////////////////////
// Sketch Retrieval Shortcuts
///////////////////////////////////////////////////////////////////////////////
vertex_sketch_raw sketch_retrieval_shortcuts::retrieve_sketch(const G &g, V v) {
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

    auto sv = retrieve_shortcuts(g, v);
    for (const auto &e : sv) {
      V tv = e.v;
      W td = d + e.d;
      auto i = pot.insert(make_pair(tv, td));
      if (!i.second && i.first->second <= td) continue;
      que.emplace(td, tv);
      i.first->second = td;
    }

    if (!FLAGS_distance_sketch_implicit_neighbor) continue;
    // TODO: faster by using only first |k| edges
    for (const auto &e : g.edges(v)) {
      W tv = to(e);
      W td = d + weight(e);
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
  return compute_sketch_retrieval_shortcuts_via_ads_unweighted(g, k, ranks, d);
}

sketch_retrieval_shortcuts compute_sketch_retrieval_shortcuts_via_ads_naive
(const G& g, size_t k, const rank_array& ranks, D d) {
  all_distances_sketches ads = compute_all_distances_sketches(g, k, ranks, d);

  sketch_retrieval_shortcuts srs;
  srs.k = ads.k;
  srs.ranks = ads.ranks;
  srs.sketches.resize(g.num_vertices());

  // (Distance, to, from)
  using tuple_type = tuple<W, V, V>;
  vector<tuple_type> es;
  {
    for (V v : g.vertices()) {
      const vertex_sketch_raw &s = ads.retrieve_sketch(g, v);
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

    vertex_sketch_raw a = srs.retrieve_sketch(g, v);
    // pretty_print(a);
    if (find(a.begin(), a.end(), entry(s, d)) != a.end()) {
//      cout << "SKIP: " << make_tuple(v, s, d) << " " << endl;
      // cout << endl;
      continue;
    }
//    cout << "IN: " << make_tuple(v, s, d) << " " << endl;
    srs.sketches[v].emplace_back(s, d);
  }

  return srs;
}

sketch_retrieval_shortcuts compute_sketch_retrieval_shortcuts_via_ads_fast
(const G& g, size_t k, const rank_array& ranks, D d) {
  all_distances_sketches ads = compute_all_distances_sketches(g, k, ranks, d);

  sketch_retrieval_shortcuts srs;
  srs.k = ads.k;
  srs.ranks = ads.ranks;
  srs.sketches.resize(g.num_vertices());

  // (Distance, to, from)
  using tuple_type = tuple<W, V, V>;
  vector<tuple_type> es;
  vector<vertex_sketch_raw> rev(g.num_vertices());
  for (V v : g.vertices()) {
    const vertex_sketch_raw &s = ads.retrieve_sketch(g, v);
    for (const auto &e : s) {
      es.emplace_back(make_tuple(e.d, e.v, v));
      rev[e.v].emplace_back(v, e.d);
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

    if (!(FLAGS_distance_sketch_implicit_neighbor && d <= 1)) {
      srs.sketches[v].emplace_back(s, d);
    }

    for (const auto &r : rev[v]) {
      // TODO: this is the bottleneck
      gs.emplace(make_tuple(r.d + d, s, r.v));  // r.v -> v -> s
    }
  }
  return srs;
}


sketch_retrieval_shortcuts compute_sketch_retrieval_shortcuts_via_ads_unweighted
(const G& g, size_t k, const rank_array& ranks, D d) {
  all_distances_sketches ads = compute_all_distances_sketches(g, k, ranks, d);

  sketch_retrieval_shortcuts srs;
  srs.k = ads.k;
  srs.ranks = ads.ranks;
  srs.sketches.resize(g.num_vertices());

  // (Distance, to, from)
  using tuple_type = tuple<W, V, V>;
  vector<tuple_type> es;
  vector<vertex_sketch_raw> rev(g.num_vertices());
  for (V v : g.vertices()) {
    const vertex_sketch_raw &s = ads.retrieve_sketch(g, v);
    for (const auto &e : s) {
      es.emplace_back(make_tuple(e.d, e.v, v));
      rev[e.v].emplace_back(v, e.d);
    }
  }
  sort(es.begin(), es.end());

  vector<vector<pair<V, V>>> gs(10);
  size_t i_es = 0;
  for (W d = 0; i_es < es.size(); ++d) {
    vector<pair<V, V>> gsd;
    if (d < (W)gs.size()) gsd = move(gs[d]);
    sort(gsd.begin(), gsd.end());
    size_t i_gsd = 0;

    for (; i_es < es.size() && get<0>(es[i_es]) == d; ++i_es) {
      const auto &e = es[i_es];
      const V s = get<1>(e);
      const V v = get<2>(e);
      if (d == 0) continue;

      while (i_gsd < gsd.size() && gsd[i_gsd] < make_pair(s, v)) ++i_gsd;
      if (i_gsd < gsd.size() && gsd[i_gsd] == make_pair(s, v)) continue;

      if (!(FLAGS_distance_sketch_implicit_neighbor && d <= 1)) {
        srs.sketches[v].emplace_back(s, d);
      }

      for (const auto &r : rev[v]) {
        W td = r.d + d;
        if (td >= (W)gs.size()) gs.resize(gs.size() * 2);
        gs[td].emplace_back(make_pair(s, r.v));  // r.v -> v -> s
      }
    }
  }
  return srs;
}

///////////////////////////////////////////////////////////////////////////////
// Update
///////////////////////////////////////////////////////////////////////////////
bool dynamic_all_distances_sketches::add_entry(V v, V s, W d) {
  size_t num_lose = 0;

  for (auto &e : ads_.sketches[v]) {
    if (e.v == s) {
      if (d >= e.d) return false;
      e.d = d;
      goto ins;
    }
    if (ads_.ranks[e.v] < ads_.ranks[s] &&
        make_pair((W)e.d, (V)e.v) < make_pair(d, s)) {
      ++num_lose;
      if (num_lose >= ads_.k) return false;
    }
  }

  assert(num_lose < ads_.k);
  ads_.sketches[v].emplace_back(s, d);

  // Remove no longer unnecessary entries
  ins:
  ads_.sketches[v] = purify_sketch(move(ads_.sketches[v]), ads_.k, ads_.ranks);
  return true;
}

void dynamic_all_distances_sketches::expand(const G &g, V v, V s, W d) {
  if (!add_entry(v, s, d)) return;

  queue<pair<V, W>> que;
  que.push({v, d});

  while (!que.empty()) {
    V x = que.front().first;
    W d = que.front().second;
    que.pop();
    for (const auto &e : g.edges(x, reverse_direction(d_))) {
      V tx = to(e);
      if (!add_entry(tx, s, d + weight(e))) continue;
      que.push({tx, d + weight(d)});
    }
  }
}

vector<V> dynamic_all_distances_sketches::shrink(const G &g, V u, V r, W durL) {
  if (u == r) return {};
  // cout << "SHRINK: " << make_tuple(u, r, durL) << endl;

  vector<V> S;
//  priority_queue<pair<W, V>, vector<pair<W, V>>, greater<pair<W, V>>> Q;
//  Q.emplace(durL, u);
  heap_.clear();
  heap_.decrease(u, durL);
  set<V> vis;  // TODO: faster
  while (!heap_.empty()) {
//    V x = Q.top().second;
//    W dxrL = Q.top().first;
//    Q.pop();
    V x = heap_.top_vertex();
    W dxrL = heap_.top_weight();
    heap_.pop();
    // cout << "DEQUE: " << make_tuple(r, x, dxrL) << endl;

    W dxrU = kInfW;
    for (const auto &e : g.edges(x, d_)) {
      W d = find_distance(ads_.sketches[to(e)], r);
      if (is_lt(d, dxrL)) dxrU = min(dxrU, d + weight(e));
    }

    if (is_lt(dxrL, dxrU)) {
      auto &sketch = ads_.sketches[x];
      {
        // cout << make_tuple(x, dxrL, r) << endl; pretty_print(sketch);
        auto ite = remove(sketch.begin(), sketch.end(), entry(r, dxrL));
        assert(ite != sketch.end());
        sketch.erase(ite, sketch.end());
      }
      S.emplace_back(x);
      for (const auto &e : g.edges(x, reverse_direction(d_))) {
        V y = to(e);
        W dyrL = find_distance(ads_.sketches[y], r);
        if (is_eq(dyrL, weight(e) + dxrL) && !vis.count(y)) {
          // Q.emplace(dyrL, y);
          heap_.decrease(y, dyrL);
          vis.insert(y);
        }
      }
    }
  }

  return S;
}

void dynamic_all_distances_sketches::re_insert(const G &g, std::vector<V> S, V r) {
  sort(S.begin(), S.end());
  heap_.clear();
  for (V x : S) {
    W dxrU = kInfW;
    for (const auto &e : g.edges(x)) {
      W dyr = find_distance(ads_.sketches[to(e)], r);
      if (dyr != kInfW) dxrU = min(dxrU, dyr + weight(e));
    }
    heap_.decrease(x, dxrU);
  }

  while (!heap_.empty()) {
    V x = heap_.top_vertex();
    W dxrU = heap_.top_weight();
    heap_.pop();

    add_entry(x, r, dxrU);  // TODO: if not, skip edge traversal?
    for (const auto &e : g.edges(x, reverse_direction(d_))) {
      V y = to(e);
      if (!binary_search(S.begin(), S.end(), y)) continue;
      W dyrU = weight(e) + dxrU;
      heap_.decrease(y, dyrU);
    }
  }
}

void dynamic_all_distances_sketches::add_edge(const G& g, V v_from, const E& e) {
  // TODO: visited flags (to avoid |add_entry|) -> lazy purify
  assert(d_ == kFwd);

  V v_to = to(e);
  for (const auto &ent : ads_.sketches[v_to]) {
    expand(g, v_from, ent.v, ent.d + weight(e));
  }
}

void dynamic_all_distances_sketches::remove_edge(const G& g, V v_from, V v_to) {
  assert(d_ == kFwd);

  vector<V> allS;
  vertex_sketch_raw sketch = ads_.sketches[v_from];
  for (const auto &ent : sketch) {
    auto S = shrink(g, v_from, ent.v, ent.d);
    allS.insert(allS.end(), S.begin(), S.end());
    re_insert(g, move(S), ent.v);
  }

  sort(allS.begin(), allS.end());
  allS.erase(unique(allS.begin(), allS.end()), allS.end());
  for (V x : allS) {
    for (const auto &edg : g.edges(x, d_)) {
      V y = to(edg);
      for (const auto &ent : ads_.sketches[y]) {
        W dxrU = weight(edg) + ent.d;
        expand(g, x, ent.v, dxrU);
      }
    }
  }
}

}  // namespace distance_sketch
}  // namespace agl
