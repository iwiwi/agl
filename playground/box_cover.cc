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
  vector<map<V, V>> X(num_v);
  vector<V> rank(num_v);
  vector<V> inv(num_v);
  vector<set<V>> A(num_v);
  for (V i = 0; i < num_v; ++i) {
    inv[i] = i;
    A[i].insert(i);
  }
  random_shuffle(inv.begin(), inv.end());
  for (int i = 0; i < num_v; ++i) {
    rank[inv[i]] = i;
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
      for (V a : A[v])
        for (V w : g.neighbors(a)) {
          // Merge & Purify
          V max_rank = X[w].rbegin()->first;
          V rv = rank[v];
          if (X[w].find(rv) != X[w].end()) continue;
          if (X[w].size() >= k && max_rank > rv) {
            X[w].erase(max_rank);
            X[w][rv] = v;
            next.insert(w);
          } else if (X[w].size() < k) {
            X[w][rv] = v;
            next.insert(w);
          }
        }
      A[v].swap(next);
    }
  }

  //
  // Select-Greedy O(n^2*k)
  //
  vector<V> centers;
  vector<bool> centered(num_v);
  map<V, V> Xs;
  priority_queue<pair<double, V>> que;

  vector<double> debug_before(num_v, 0.0);
  while (estimated_cardinality(g, Xs, k) < max(num_v, k - 1)) {
    V selected_v = -1;

    double ec_before = estimated_cardinality(g, Xs, k);
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

      double ec_delta = ec_tmp - ec_before;
      if (debug_before[v] > 0.0 && debug_before[v] < ec_delta) {
        cerr << v << endl;
        cerr << debug_before[v] << " " << ec_delta << endl;
        for (auto p : X[v]) {
          cerr << "X " << p.first << " " << p.second << endl;
        }
        int d = 0;
        V dd = 0;
        for (auto p : Xs) {
          d++;
          if (d == k) dd = p.first;
          cerr << "Xs" << p.first << " " << p.second << endl;
        }
        int c = 0;
        V cc = 0;
        for (auto p : tmp) {
          c++;
          if (c == k) cc = p.first;
          cerr << "T " << p.first << " " << p.second << endl;
        }
        cerr << dd << " " << cc << endl;
        assert(false);
      }

      debug_before[v] = ec_delta;
      if (v == 11) {
        for (auto p : Xs) {
          cerr << "Xs " << p.first << " " << p.second << endl;
        }
        for (auto p : tmp) {
          cerr << "tmp " << p.first << " " << p.second << endl;
        }
      }
    }
    cerr << "Naive " << selected_v << " " << argmax << endl;

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