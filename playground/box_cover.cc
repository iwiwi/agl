#include "box_cover.h"
#include <sys/time.h>

vector<V> box_cover_memb(const G &g, W radius) {
  V num_v = g.num_vertices();
  vector<vector<pair<V, W>>> node_lists;
  map<size_t, set<V>> excluded_mass_map;
  {
    vector<pair<size_t, V>> center_candidates;
    for (V pv = 0; pv < num_v; ++pv) {
      queue<pair<V, W>> que;
      vector<bool> vis(num_v, false);
      vector<pair<V, W>> nodes;

      que.push(make_pair(pv, 0));
      nodes.push_back(make_pair(pv, 0));
      vis[pv] = true;

      while (!que.empty()) {
        V v = que.front().first;
        W dist = que.front().second;
        que.pop();
        if (dist >= radius) continue;
        for (V u : g.neighbors(v)) {
          if (vis[u]) continue;
          que.push(make_pair(u, dist + 1));
          nodes.push_back(make_pair(u, dist + 1));
          vis[u] = true;
        }
      }

      node_lists.push_back(nodes);
      center_candidates.push_back(make_pair(nodes.size(), pv));
    }
    for (pair<size_t, V> p : center_candidates) {
      excluded_mass_map[p.first].insert(p.second);
    }
  }

  set<V> covered_nodes;
  set<V> center_nodes;
  vector<W> central_distance(num_v, num_v);
  while (covered_nodes.size() < (size_t)num_v) {
    V center_node_found;
    while (true) {
      V node, maximum_key;
      while (true) {
        maximum_key = excluded_mass_map.rbegin()->first;
        set<V> &nodes = excluded_mass_map.rbegin()->second;
        auto it = nodes.begin();
        advance(it, agl::random(nodes.size()));
        node = *it;
        if (center_nodes.find(node) != center_nodes.end()) {
          nodes.erase(it);
          if (nodes.empty()) excluded_mass_map.erase(maximum_key);
        } else {
          break;
        }
      }

      V mass = 0;
      for (pair<V, W> p : node_lists[node]) {
        if (covered_nodes.find(p.first) == covered_nodes.end()) {
          mass++;
        }
      }
      excluded_mass_map[maximum_key].erase(node);
      if (excluded_mass_map[maximum_key].empty())
        excluded_mass_map.erase(maximum_key);
      if (mass == maximum_key) {
        center_node_found = node;
        break;
      } else {
        excluded_mass_map[mass].insert(node);
      }
    }
    center_nodes.insert(center_node_found);
    for (pair<V, W> p : node_lists[center_node_found]) {
      V i = p.first;
      W d = p.second;
      covered_nodes.insert(i);
      central_distance[i] = max(central_distance[i], d);
    }
  }

  vector<int> ret(center_nodes.begin(), center_nodes.end());
  return ret;
}

vector<V> box_cover_burning(const G &g, W radius) { return {}; }

double GetCurrentTimeSec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

vector<map<V, V>> naive_build_sketch(const G &g, const W radius, const int k,
                                     const vector<V> &rank) {
  V num_v = g.num_vertices();
  vector<map<V, V>> naive_X(num_v);
  for (V i = 0; i < num_v; ++i) {
    map<V, V> tmp;
    vector<bool> vis(num_v, false);
    queue<pair<V, W>> que;
    que.push(make_pair(i, 0));
    vis[i] = true;
    tmp[rank[i]] = i;

    while (!que.empty()) {
      V v = que.front().first;
      W dist = que.front().second;
      que.pop();
      if (dist == radius) break;
      for (V u : g.neighbors(v)) {
        if (vis[u]) continue;
        que.push(make_pair(u, dist + 1));
        vis[u] = true;
        tmp[rank[u]] = u;
      }
    }
    int cnt = 0;
    for (pair<V, V> p : tmp) {
      naive_X[i][p.first] = p.second;
      cnt++;
      if (cnt == k) break;
    }
  }

  return naive_X;
}

vector<map<V, V>> build_sketch(const G &g, const W radius, const int k,
                               const vector<V> &rank, const vector<V> &inv) {
  V num_v = g.num_vertices();
  vector<map<V, V>> X(num_v);
  vector<set<V>> previous_added(num_v);
  for (V i = 0; i < num_v; ++i) {
    previous_added[i].insert(i);
  }
  //
  // Build-Sketches O((n+m)*rad)
  //
  for (V i = 0; i < num_v; ++i) {
    X[i][rank[i]] = i;
  }
  for (W d = 0; d < radius; ++d) {
    for (V v : inv) {
      set<V> next;
      for (V a : previous_added[v])
        for (V neighbor : g.neighbors(a)) {
          // Merge & Purify
          V max_rank = X[neighbor].rbegin()->first;
          V rv = rank[v];
          if (X[neighbor].find(rv) != X[neighbor].end()) continue;
          if (X[neighbor].size() >= k && max_rank > rv) {
            X[neighbor].erase(max_rank);
            X[neighbor][rv] = v;
            next.insert(neighbor);
          } else if (X[neighbor].size() < k) {
            X[neighbor][rv] = v;
            next.insert(neighbor);
          }
        }

      previous_added[v].swap(next);
    }
  }
  return X;
}

double estimated_cardinality(const G &g, const map<V, V> X, const int k) {
  if (X.size() < k) {
    return X.size();
  }
  int cnt = 0;
  for (pair<V, V> p : X) {
    cnt++;
    if (cnt == k) {
      return (double)(k - 1) / p.first * g.num_vertices();
    }
  }
  assert(false);
  return (k - 1);
}

vector<V> box_cover_sketch(const G &g, W radius) {
  const int k = 200;
  const V num_v = g.num_vertices();
  vector<V> rank(num_v);
  vector<V> inv(num_v);
  for (V i = 0; i < num_v; ++i) {
    inv[i] = i;
  }
  random_shuffle(inv.begin(), inv.end());
  for (int i = 0; i < num_v; ++i) {
    rank[inv[i]] = i;
  }

  //
  // Build-Sketches O((n+m)*rad)
  //
  vector<map<V, V>> X = build_sketch(g, radius, k, rank, inv);

  //
  // Select-Greedy O(n^2*k)
  //
  vector<V> centers;
  vector<bool> centered(num_v);
  map<V, V> Xs;

  vector<double> debug_before(num_v, 0.0);
  while (estimated_cardinality(g, Xs, k) < max(num_v, k - 1)) {
    V selected_v = -1;

    double argmax = 0.0;
    for (V v = 0; v < num_v; v++) {
      if (centered[v]) continue;
      map<V, V> tmp(Xs);
      tmp.insert(X[v].begin(), X[v].end());

      double ec_tmp = estimated_cardinality(g, tmp, k);
      if (argmax < ec_tmp) {
        argmax = ec_tmp;
        selected_v = v;
      }
    }

    assert(selected_v >= 0);

    centers.push_back(selected_v);
    centered[selected_v] = true;
    if (Xs.size() == 0) {
      for (pair<V, V> p : X[selected_v]) {
        Xs[p.first] = p.second;
        if (Xs.size() == k) goto purified;
      }
    }
    // Merge-and-Purify
    for (pair<V, V> p : X[selected_v]) {
      V max_rank = Xs.rbegin()->first;
      if (Xs.size() > k && p.first > max_rank) break;
      if (Xs.size() > k) {
        Xs.erase(max_rank);
      }
      Xs[p.first] = p.second;
    }
  purified:
    ;
  }

  return centers;
}