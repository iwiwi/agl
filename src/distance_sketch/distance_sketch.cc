#include "distance_sketch.h"
#include "agl.h"
using namespace std;

DEFINE_int32(distance_sketch_k, 16, "");
DEFINE_bool(distance_sketch_implicit_neighbor, false, "");
DEFINE_int32(distance_sketch_srs_cache_size, 15, "");

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


///////////////////////////////////////////////////////entry////////////////////////
// Static Sketch
///////////////////////////////////////////////////////////////////////////////
vertex_sketch_raw sort_by_vertices(vertex_sketch_raw sketch) {
  sort(sketch.begin(), sketch.end(),
       [](entry e1, entry e2) { return e1.v < e2.v; });
  return sketch;
}

vertex_sketch_raw sort_by_distances(vertex_sketch_raw sketch) {
  sort(sketch.begin(), sketch.end(),
       [](entry e1, entry e2) {
    return e1.d != e2.d ? e1.d < e2.d : e1.v < e2.v;
  });
  return sketch;
}

vertex_sketch_raw sort_by_ranks(vertex_sketch_raw sketch, const rank_array &ranks) {
  sort(sketch.begin(), sketch.end(),
       [&](entry e1, entry e2) { return ranks[e1.v] < ranks[e2.v]; });
  return sketch;
}

vertex_sketch_raw purify_sketch(vertex_sketch_raw sketch, size_t k, const rank_array &ranks) {
  sketch = sort_by_distances(move(sketch));

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

  return sort_by_vertices(move(new_sketch));
}

W find_distance(const vertex_sketch_raw &sketch, V v) {
  if (sketch.empty()) return kInfW;
  size_t l = 0, r = sketch.size();
  while (r - l > 1) {
    size_t m = (r + l) / 2;
    if (v < sketch[m].v) r = m;
    else l = m;
  }
  if (sketch[l].v == v) return sketch[l].d;
  else return kInfW;
}

vertex_sketch_raw remove_entry(vertex_sketch_raw sketch, V v) {
  size_t i = 0;
  for (; i < sketch.size(); ++i) if (sketch[i].v == v) break;
  if (i == sketch.size()) return sketch;
  for (; i + 1 < sketch.size(); ++i) sketch[i] = sketch[i + 1];
  sketch.pop_back();
  return move(sketch);
}

vertex_sketch_raw insert_entry(vertex_sketch_raw sketch, V v, W d) {
  auto i = sketch.begin();
  for (; i != sketch.end(); ++i) if (i->v >= v) break;
  if (i != sketch.end() && i->v == v) {
    i->d = d;
  } else {
    sketch.insert(i, entry(v, d));
  }
  return move(sketch);
}

vertex_sketch_raw compute_all_distances_sketch_from
(const G& g, V v, size_t k, const std::vector<rank_type> &ranks, D d) {
  vertex_sketch_raw sketch;

  auto f = [&](V v, W w) -> bool {
    sketch.emplace_back(v, w);
    return true;
  };
  visit_by_distance(g, v, f, d);

  return purify_sketch(move(sketch), k, ranks);
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
  for (size_t i : make_irange(ads.sketches.size())) {
    ofs << i << ":";
    pretty_print(ads.sketches[i], ofs);
  }
}


///////////////////////////////////////////////////////////////////////////////
// Estimation
///////////////////////////////////////////////////////////////////////////////
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
  if (pot_.size() < (size_t)g.num_vertices()) {
    pot_.resize(max((pot_.size() + 10) * 2, (size_t)g.num_vertices()), kInfW);
  }

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

  priority_queue<entry, vector<entry>, entry_comparator_greater_distance> que;
  vector<V> vis;
  auto enqueue = [&](V tv, W td) -> void {
    if (td >= pot_[tv]) return;
    que.emplace(tv, td);
    if (pot_[tv] == kInfW) vis.emplace_back(tv);
    pot_[tv] = td;
  };

  enqueue(v, 0);
  while (!que.empty()) {
    V v = que.top().v;
    W d = que.top().d;
    que.pop();
    if (d > pot_[v] || !add_entry(v, d)) continue;

    auto sv = retrieve_shortcuts(g, v);
    for (const auto &e : sv) {
      enqueue(e.v, d + e.d);
    }

    if (!FLAGS_distance_sketch_implicit_neighbor) continue;
    for (const auto &e : g.edges(v)) {
      enqueue(to(e), d + weight(e));
    }
  }

  for (V v : vis) pot_[v] = kInfW;
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

  for (auto &s : srs.sketches) sort(s.begin(), s.end());
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


  for (auto &s : srs.sketches) sort(s.begin(), s.end());
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
  for (auto &s : srs.sketches) sort(s.begin(), s.end());
  return srs;
}

///////////////////////////////////////////////////////////////////////////////
// ADS Update
///////////////////////////////////////////////////////////////////////////////

bool dynamic_all_distances_sketches::add_entry(V v, V s, W d) {
  size_t num_lose = 0;
  auto &a = ads_.sketches[v];
  if (find_distance(a, s) <= d) return false;

  const size_t n = a.size();
  size_t i_ins = n;
  for (size_t i = 0; i < n; ++i) {
    auto &e = a[i];

    if (e.v == s) {
      if (d >= e.d) {
        pretty_print(a);
        assert(false);  // return false;
      }
      e.d = d;
      goto ins;
    } else if (e.v > s && i_ins == n) {
      i_ins = i;
    }

    if (ads_.ranks[e.v] < ads_.ranks[s] &&
        make_pair((W)e.d, (V)e.v) < make_pair(d, s)) {
      ++num_lose;
      if (num_lose >= ads_.k) return false;
    }
  }

  assert(num_lose < ads_.k);
  a.insert(a.begin() + i_ins, entry(s, d));

  ins:
  vs_dirty_ads_.emplace_back(v);  // Lazy purification
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

  vector<V> S;
  heap_.clear();
  heap_.decrease(u, durL);
  set<V> vis;  // TODO: faster
  while (!heap_.empty()) {
    V x = heap_.top_vertex();
    W dxrL = heap_.top_weight();
    heap_.pop();

    W dxrU = kInfW;
    for (const auto &e : g.edges(x, d_)) {
      W d = find_distance(ads_.sketches[to(e)], r);
      if (is_lt(d, dxrL)) dxrU = min(dxrU, d + weight(e));
    }

    if (is_lt(dxrL, dxrU)) {
      auto &sketch = ads_.sketches[x];
      {
        auto ite = remove(sketch.begin(), sketch.end(), entry(r, dxrL));
        assert(ite != sketch.end());
        sketch.erase(ite, sketch.end());
      }
      S.emplace_back(x);
      for (const auto &e : g.edges(x, reverse_direction(d_))) {
        V y = to(e);
        W dyrL = find_distance(ads_.sketches[y], r);
        if (is_eq(dyrL, weight(e) + dxrL) && !vis.count(y)) {
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

    if (!add_entry(x, r, dxrU)) continue;
    for (const auto &e : g.edges(x, reverse_direction(d_))) {
      V y = to(e);
      if (!binary_search(S.begin(), S.end(), y)) continue;
      W dyrU = weight(e) + dxrU;
      heap_.decrease(y, dyrU);
    }
  }
}

void dynamic_all_distances_sketches::add_edge(const G& g, V v_from, const E& e) {
  assert(d_ == kFwd);
  assert(vs_dirty_ads_.empty());

  V v_to = to(e);
  for (const auto &ent : ads_.sketches[v_to]) {
    expand(g, v_from, ent.v, ent.d + weight(e));
  }

  purify_lazy();
}

void dynamic_all_distances_sketches::remove_edge(const G& g, V v_from, V v_to) {
  assert(d_ == kFwd);
  assert(vs_dirty_ads_.empty());

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

  purify_lazy();
}

void dynamic_all_distances_sketches::purify_lazy() {
  sort(vs_dirty_ads_.begin(), vs_dirty_ads_.end());
  vs_dirty_ads_.erase(unique(vs_dirty_ads_.begin(), vs_dirty_ads_.end()), vs_dirty_ads_.end());
  for (V v : vs_dirty_ads_) {
    auto &a = ads_.sketches[v];
    a = purify_sketch(move(a), ads_.k, ads_.ranks);
  }
  vs_dirty_ads_.clear();
}

///////////////////////////////////////////////////////////////////////////////
// SRS Update
///////////////////////////////////////////////////////////////////////////////
bool dynamic_sketch_retrieval_shortcuts::add_entry(const G &g, V v, V s, W d) {
  size_t num_lose = 0;
  load_cache(g, v);
  auto &a = ads_caches_[v].ads;

  for (auto &e : a) {
    if (e.v == s) {
      if (d >= e.d) return false;
      e.d = d;
      goto ins;
    }
    if (srs_.ranks[e.v] < srs_.ranks[s] &&
        make_pair((W)e.d, (V)e.v) < make_pair(d, s)) {
      ++num_lose;
      if (num_lose >= srs_.k) return false;
    }
  }

  assert(num_lose < srs_.k);
  a.emplace_back(s, d);

  // Remove no longer unnecessary entries
  ins:
  a = purify_sketch(move(a), srs_.k, srs_.ranks);
  ads_caches_[v].is_dirty = true;
  return true;
}

void dynamic_sketch_retrieval_shortcuts::expand(const G &g, V v, V s, W d) {
  load_cache(g, v);
  if (!add_entry(g, v, s, d)) return;

  queue<pair<V, W>> que;
  que.push({v, d});

  while (!que.empty()) {
    const V x = que.front().first;
    const W d = que.front().second;
    que.pop();

    // cout << "EXPAND: "<< make_tuple(d, s, x) << endl;
    {
      auto i = srs_invalidation_.insert(make_pair(x, make_pair(d, s)));
      if (!i.second) i.first->second = min(i.first->second, make_pair(d, s));
    }

    for (const auto &e : g.edges(x, reverse_direction(d_))) {
      const V tx = to(e);
      const W td = d + weight(e);
      if (!add_entry(g, tx, s, d + weight(e))) continue;
      que.push({tx, td});
    }
  }
}



void dynamic_sketch_retrieval_shortcuts::add_edge(const G& g, V v_from, const E& e) {
  // TODO: visited flags (to avoid |add_entry|) -> lazy purify
  assert(d_ == kFwd);

  const V v_to = to(e);
  cout << v_from << "===>" << v_to << endl;
  load_cache(g, v_to);
  JLOG_ADD_BENCHMARK("expand") {
  for (const auto &ent : ads_caches_[v_to].ads) {
    expand(g, v_from, ent.v, ent.d + weight(e));
  }
  }

  JLOG_ADD_BENCHMARK("srs") {
  assert(srs_new_sketches_.empty());
  for (const auto &i : srs_invalidation_) {
    const V v = i.first;
    // Removing invalidated SRS entries
    {
      auto ite = srs_new_sketches_.emplace(v, srs_.sketches[v]);
      auto &s = ite.first->second;
      auto si = remove_if(s.begin(), s.end(), [&](const entry &e) {
        return make_pair(e.d, e.v) >= i.second;
      });
      s.erase(si, s.end());
    }
    // Pushing ADS entries to be examined
    for (const auto &e : ads_caches_[v].ads) {
      if (make_pair(e.d, e.v) >= i.second) {
        srs_tentative_entry_que_.emplace(e.d, e.v, v);
      }
    }
  }
  srs_invalidation_.clear();

  srs_propagation_distance_.clear();
  while (!srs_tentative_entry_que_.empty()) {
    auto ue = srs_tentative_entry_que_.top();
    srs_tentative_entry_que_.pop();

    if (!srs_tentative_entry_que_.empty() && srs_tentative_entry_que_.top() == ue) {
      srs_tentative_entry_que_.pop();
    }

    const W d = get<0>(ue);
    const V s = get<1>(ue);
    const V v = get<2>(ue);

    // (s, d) is in SRS[v] now?
    W prv_d = find_distance(srs_.sketches[v], s);

    // Check if (s, d) is necessary for SRS[v]
    W new_d = d;
    load_cache(g, v);
    const auto &a = ads_caches_[v].ads;
    for (const auto &ae : a) {
      if (ae.v == v || ae.v == s) continue;
      const V tv = ae.v;
      const W td = ae.d;
      if (is_le(d, td)) continue;

      //const W td2 = find_distance(srs_.sketches[tv], s);
      W td2;
      if (srs_new_sketches_.count(tv)) td2 = find_distance(srs_new_sketches_[tv], s);
      else td2 = find_distance(srs_.sketches[tv], s);

      if (td2 == kInfW || td + td2 != d) continue;
      assert(td + td2 >= d);

      new_d = kInfW;  // not necessary
      break;
    }

    if (prv_d != new_d) {
      // cout << "CHANGE!: " << v << " -> " << s << ": " << prv_d << "->" << new_d << endl;
      // TODO: Need to propagate
      propagate(g, v, s, d);
    }

    if (new_d != kInfW) {
      // Retrieval path not found
      assert(srs_new_sketches_.count(v));  // TODO: remove this line
      auto &srs = srs_new_sketches_[v];
      srs = insert_entry(move(srs), s, new_d);
    }
  }

  // Finalize
  for (auto &i : srs_new_sketches_) srs_.sketches[i.first] = move(i.second);
  srs_new_sketches_.clear();
  cout << num_cached_vertices_ << " / " << g.num_vertices() << endl;
  purify_cache(g, 0);

  }
}

void dynamic_sketch_retrieval_shortcuts::propagate(const G& g, V v, V s, W d) {
  queue<pair<V, W>> que;

  auto enque = [&](V tv, W td) -> void {
    if (srs_propagation_distance_.count({s, tv})) return;

    load_cache(g, tv);
    if (find_distance(ads_caches_[tv].ads, s) != td) return;

    que.push({tv, td});
    srs_propagation_distance_.insert(make_pair(make_pair(s, tv), td));
  };

  for (const auto &e : g.edges(v, reverse_direction(d_))) {
    enque(to(e), d + weight(e));
  }

  while (!que.empty()) {
    const V x = que.front().first;
    const W d = que.front().second;
    que.pop();

    // Invalidation
    {
      auto i = srs_new_sketches_.insert(make_pair(x, vertex_sketch_raw()));
      auto &srs = i.first->second;
      if (i.second) srs = srs_.sketches[x];
      srs = remove_entry(move(srs), s);
      srs_tentative_entry_que_.emplace(d, s, x);
    }

    // Edge traversal
    for (const auto &e : g.edges(x, reverse_direction(d_))) {
      enque(to(e), d + weight(e));
    }
  }
}
}  // namespace distance_sketch
}  // namespace agl

