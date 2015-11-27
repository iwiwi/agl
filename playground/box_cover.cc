#include "box_cover.h"

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

vector<V> box_cover_sketch(const G &g, W radius, const int k) {
  V num_v = g.num_vertices();
  vector<V> rank(num_v);
  vector<map<V, V>> X(num_v);
  vector<V> inv(num_v);
  vector<queue<V>> A(num_v);
  for (V i = 0; i < num_v; ++i) {
    inv[i] = i;
    A[i].push(i);
  }
  random_shuffle(inv.begin(), inv.end());
  for (int i = 0; i < num_v; ++i) {
    rank[inv[i]] = i;
  }
  //
  // Build-Sketches
  //
  for (V i = 0; i < num_v; ++i) {
    X[i][rank[i]] = i;
  }

  for (W d = 0; d < radius; ++d) {
    for (V v : inv) {
      size_t size_A = A[v].size();
      for (int q = 0; q < size_A; ++q) {
        V a = A[v].front();
        A[v].pop();
        for (V w : g.neighbors(a)) {
          // Merge & Purify
          V max_rank = X[w].rbegin()->first;
          V rv = rank[v];
          if (X[w].size() >= k && max_rank > rv && X[w].find(rv) == X[w].end()) {
            X[w].erase(max_rank);
            X[w][rv] = v;
            A[v].push(w);
          } else if (X[w].size() < k && X[w].find(rv) == X[w].end()) {
            X[w][rv] = v;
            A[v].push(w);
          }
        }
      }
    }
  }

  //
  // Select-Greedy O(n^2*k)
  //
  vector<V> centers;
  map<V, V> Xs;

  return {};
};