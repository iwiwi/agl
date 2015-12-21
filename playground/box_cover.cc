#include "box_cover.h"











double coverage(const G &g, const vector<V> &s, W rad,
                vector<bool> &is_covered) {
  vector<V> covered;
  vector<bool> vis(g.num_vertices());
  queue<pair<V, W>> que;
  for (V center : s) {
    que.push(make_pair(center, 0));
    covered.push_back(center);
    vis[center] = true;
  }
  while (!que.empty()) {
    V v = que.front().first;
    W dist = que.front().second;
    que.pop();
    if (dist == rad) break;
    for (V u : g.neighbors(v)) {
      if (vis[u]) continue;
      vis[u] = true;
      que.push(make_pair(u, dist + 1));
      covered.push_back(u);
    }
  }
  assert(covered.size() <= (size_t)g.num_vertices());

  for (V c : covered) {
    is_covered[c] = true;
  }

  double ret = covered.size();
  return ret / g.num_vertices();
}

double coverage(const G &g, const vector<V> &s, W rad) {
  vector<bool> dummy(g.num_vertices(), false);
  return coverage(g, s, rad, dummy);
}

vector<V> merge_and_purify(set<V> &parent, const vector<V> &sorted_vec,
                           const int k) {
  vector<V> delta;
  for (V p_rank : sorted_vec) {
    if (parent.size() == (size_t)k && p_rank > *parent.rbegin()) break;
    size_t prev = parent.size();
    parent.insert(p_rank);
    if (parent.size() > prev) delta.push_back(p_rank);
    while (parent.size() > (size_t)k) parent.erase(*parent.rbegin());
  }
  return delta;
}

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

  vector<V> ret(center_nodes.begin(), center_nodes.end());
  return ret;
}

void burning_splitted(const G &g, W radius, set<V> &solution,
                      vector<vector<V>> &boxes) {
  V num_v = g.num_vertices();
  V prev_size = -1;
  while (prev_size < (V)solution.size()) {
    prev_size = solution.size();
    //
    // Remove unnecessary boxes.
    //
    {
      for (int i = 0; i < num_v; ++i) {
        sort(boxes[i].begin(), boxes[i].end());
        boxes[i].erase(unique(boxes[i].begin(), boxes[i].end()),
                       boxes[i].end());
      }

      for (V i = 0; i < num_v; ++i)
        for (V j = i + 1; j < num_v; ++j) {
          size_t cnt = 0;
          for (size_t ci = 0, cj = 0;
               ci < boxes[i].size() && cj < boxes[j].size();) {
            if (boxes[i][ci] == boxes[j][cj]) {
              cnt++;
              ci++;
              cj++;
            } else if (boxes[i][ci] < boxes[j][cj]) {
              ci++;
            } else {
              cj++;
            }
          }
          if (cnt == boxes[j].size()) {
            boxes[j].clear();
          } else if (cnt == boxes[i].size()) {
            boxes[i].clear();
          }
        }
    }

    //
    // Remove unnecessary nodes.
    //
    {
      vector<vector<V>> containd_box_list(num_v);
      for (V i = 0; i < num_v; ++i)
        for (V v : boxes[i]) containd_box_list[v].push_back(i);

      for (V i = 0; i < num_v; ++i)
        sort(containd_box_list[i].begin(), containd_box_list[i].end());

      for (V i = 0; i < num_v; ++i)
        for (V j = i + 1; j < num_v; ++j) {
          if (containd_box_list[i].size() == 0 ||
              containd_box_list[j].size() == 0)
            continue;
          size_t cnt = 0;
          for (size_t ci = 0, cj = 0; ci < containd_box_list[i].size() &&
                                          cj < containd_box_list[j].size();) {
            if (containd_box_list[i][ci] == containd_box_list[j][cj]) {
              cnt++;
              ci++;
              cj++;
            } else if (containd_box_list[i][ci] < containd_box_list[j][cj]) {
              ci++;
            } else {
              cj++;
            }
          }

          if (cnt == containd_box_list[i].size()) {
            containd_box_list[j].clear();
          } else if (cnt == containd_box_list[j].size()) {
            containd_box_list[i].clear();
          }
        }

      boxes.clear();
      boxes.resize(num_v);
      for (V v = 0; v < num_v; ++v)
        for (V b : containd_box_list[v]) boxes[b].push_back(v);

      // Remove pairs of unnecessary twin boxes
      for (V i = 0; i < num_v; ++i)
        for (V j = i + 1; j < num_v; ++j)
          if (containd_box_list[i].size() == 2 &&
              containd_box_list[j].size() == 2) {
            vector<V> &bi1 = boxes[containd_box_list[i][0]];
            vector<V> &bi2 = boxes[containd_box_list[i][1]];
            vector<V> &bj1 = boxes[containd_box_list[j][0]];
            vector<V> &bj2 = boxes[containd_box_list[j][1]];
            if (bi1.size() != 2 || bi2.size() != 2 || bj1.size() != 2 ||
                bj2.size() != 2)
              continue;

            V k1 = bi1[0] == i ? bi1[1] : bi1[0];
            V k2 = bi2[0] == i ? bi2[1] : bi2[0];
            V l1 = bj1[0] == j ? bj1[1] : bj1[0];
            V l2 = bj2[0] == j ? bj2[1] : bj2[0];

            if (k1 == l1 && k2 == l2) {
              bi2.clear();
              bj1.clear();
            } else if (k1 == l2 && k2 == l1) {
              bi2.clear();
              bj2.clear();
            }
          }
    }

    //
    // Search for boxes that must be contained in the solution.
    //
    {
      vector<vector<V>> containd_box_list(num_v);
      for (V i = 0; i < num_v; ++i)
        for (V v : boxes[i]) containd_box_list[v].push_back(i);
      for (V v = 0; v < num_v; ++v)
        if (containd_box_list[v].size() == 1) {
          V b = containd_box_list[v][0];
          solution.insert(b);
          vector<V> &covered = boxes[b];
          for (V c : covered) containd_box_list[c].clear();
        }

      boxes.clear();
      boxes.resize(num_v);
      for (V v = 0; v < num_v; ++v)
        for (V b : containd_box_list[v]) boxes[b].push_back(v);
    }
  }
}

vector<V> box_cover_burning(const G &g, W radius) {
  if (radius == 0) {
    vector<V> ret;
    for (V v = 0; v < g.num_vertices(); v++) ret.push_back(v);
    return ret;
  }
  set<V> solution;

  //
  // Create all possible boxes.
  //
  V num_v = g.num_vertices();
  vector<vector<V>> boxes(num_v);
  {
    for (V i = 0; i < num_v; ++i) {
      vector<bool> vis(num_v);
      queue<pair<V, W>> que;

      que.push(make_pair(i, 0));
      vis[i] = true;
      boxes[i].push_back(i);

      while (!que.empty()) {
        V q = que.front().first;
        W dist = que.front().second;
        que.pop();
        if (dist == radius) break;
        for (V v : g.neighbors(q)) {
          if (vis[v]) continue;
          que.push(make_pair(v, dist + 1));
          vis[v] = true;
          boxes[i].push_back(v);
        }
      }
      sort(boxes[i].begin(), boxes[i].end());
    }
  }

  //
  // System split.
  // Find the best solution by DFS
  // The box-covering for tree networks could be performed in O(N^3)
  // while for regular networks it requires O(2^N).
  //
  deque<pair<vector<vector<V>>, set<V>>> que;
  que.push_front(make_pair(boxes, solution));
  size_t min_size = num_v;
  while (!que.empty()) {
    vector<vector<V>> qboxes = que.front().first;
    set<V> qsolution = que.front().second;
    que.pop_front();
    if (qsolution.size() >= min_size) continue;

    burning_splitted(g, radius, qsolution, qboxes);

    vector<V> s(qsolution.begin(), qsolution.end());
    double qcoverage = coverage(g, s, radius);
    if (qcoverage == 1.0) {
      if (min_size > qsolution.size()) {
        solution = qsolution;
        min_size = solution.size();
      }
      continue;
    }

    // Find the node that is in the smallest number of boxes
    vector<set<V>> containd_box_list(num_v);
    vector<V> covered_largest(num_v, 0);
    for (V i = 0; i < num_v; ++i)
      for (V v : qboxes[i]) {
        containd_box_list[v].insert(i);
        if ((size_t)covered_largest[v] < qboxes[i].size())
          covered_largest[v] = qboxes[i].size();
      }
    V selected_v = -1;
    size_t min_list_size = num_v;
    V max_covered = 0;
    for (int v = 0; v < num_v; ++v) {
      if (containd_box_list[v].empty()) continue;
      if (containd_box_list[v].size() < min_list_size) {
        min_list_size = containd_box_list[v].size();
        max_covered = covered_largest[v];
        selected_v = v;
      } else if (containd_box_list[v].size() == min_list_size &&
                 covered_largest[v] > max_covered) {
        min_list_size = containd_box_list[v].size();
        max_covered = covered_largest[v];
        selected_v = v;
      }
    }

    assert(selected_v >= 0 && containd_box_list[selected_v].size() >= 2);

    // After Selected
    for (V onlycontain : containd_box_list[selected_v]) {
      // Generate subboxes
      vector<vector<V>> subboxes(num_v);
      for (V v = 0; v < num_v; v++) {
        if (v != onlycontain &&
            containd_box_list[selected_v].find(v) !=
                containd_box_list[selected_v].end())
          continue;
        subboxes[v].assign(qboxes[v].begin(), qboxes[v].end());
      }
      que.push_front(make_pair(subboxes, qsolution));
    }
  }
  vector<V> ret(solution.begin(), solution.end());
  return ret;
}

// Naive BFS method of Build-Sketch
vector<vector<V>> naive_build_sketch(const G &g, const W radius, const int k,
                                     const vector<V> &rank,
                                     const vector<V> &inv,
                                     const vector<bool> &is_covered) {
  V num_v = g.num_vertices();
  vector<vector<V>> naive_X(num_v);
  for (V i = 0; i < num_v; ++i) {
    set<V> tmp;
    vector<bool> vis(num_v, false);
    queue<pair<V, W>> que;
    que.push(make_pair(i, 0));
    vis[i] = true;
    tmp.insert(rank[i]);

    while (!que.empty()) {
      V v = que.front().first;
      W dist = que.front().second;
      que.pop();
      if (dist == radius) break;
      for (V u : g.neighbors(v)) {
        if (vis[u]) continue;
        que.push(make_pair(u, dist + 1));
        vis[u] = true;
        tmp.insert(rank[u]);
      }
    }
    int cnt = 0;
    for (V p : tmp) {
      if (is_covered[inv[p]]) continue;
      naive_X[i].push_back(p);
      cnt++;
      if (cnt == k) break;
    }
  }

  return naive_X;
}

vector<vector<V>> build_sketch(const G &g, const W radius, const int k,
                               const vector<V> &rank, const vector<V> &inv,
                               const vector<bool> &is_covered) {
  V num_v = g.num_vertices();
  vector<set<V>> X(num_v);
  vector<set<V>> previous_added(num_v);
  //
  // Build-Sketches O((n+m)*rad)
  //
  for (V i = 0; i < num_v; ++i) {
    if (is_covered[i]) continue;
    X[i].insert(rank[i]);
    previous_added[i].insert(i);
  }

  for (W d = 0; d < radius; ++d) {
    for (V v : inv) {
      set<V> next;
      for (V a : previous_added[v])
        for (V neighbor : g.neighbors(a)) {
          // Merge & Purify
          V rv = rank[v];
          if (X[neighbor].find(rv) != X[neighbor].end()) continue;
          if (X[neighbor].size() >= (size_t)k) {
            V max_rank = *X[neighbor].rbegin();
            if (max_rank <= rv) continue;
            X[neighbor].erase(max_rank);
            X[neighbor].insert(rv);
            next.insert(neighbor);
          } else {
            X[neighbor].insert(rv);
            next.insert(neighbor);
          }
        }

      previous_added[v].swap(next);
    }
  }

  vector<vector<V>> ret;
  for (V v = 0; v < num_v; ++v) {
    vector<V> sketch(X[v].begin(), X[v].end());
    ret.push_back(sketch);
  }

  return ret;
}

void select_greedily(const G &g, const vector<vector<V>> &X, vector<V> &centers,
                     vector<bool> &centered, const int k) {
  assert(g.num_vertices() > k);
  //
  // Variables
  //
  V num_v = g.num_vertices();
  set<V> Xs;
  priority_queue<pair<V, V>, vector<pair<V, V>>, greater<pair<V, V>>> que[2];
  vector<multimap<V, V>> T(k + 2);
  vector<V> k1 = vector<V>(num_v);
  vector<V> k2 = vector<V>(num_v);
  vector<V> c = vector<V>(num_v);
  vector<V> last_blue(num_v);
  vector<int> type(num_v, 0);
  vector<set<V>> I(num_v);
  vector<bool> removed(num_v, false);
  vector<bool> covered_rank(num_v, false);

  auto insert_as_type = [&](V box, int target_type) {
    if (target_type == 0) {
      type[box] = 0;
      T[k2[box] + 1].insert({last_blue[box], box});
      que[0].push({last_blue[box], box});
    } else {
      type[box] = 1;
      T[k2[box]].insert({last_blue[box], box});
      que[1].push({k2[box], box});
    }
  };
  auto remove_covered_ranks = [&](V box) {
    while (covered_rank[last_blue[box]]) {
      k1[box]--;
      c[box]--;
      if (k1[box] == 0) {
        removed[box] = true;
        break;
      }
      last_blue[box] = X[box][k1[box] - 1];
    }
  };
  auto remove_multimap_pair = [&](V from, const pair<V, V> &p) {
    auto it_p = T[from].equal_range(p.first);
    for (auto it = it_p.first; it != T[from].end(); it++) {
      if ((*it).second == p.second) {
        T[from].erase(it);
        break;
      }
      if (it == it_p.second) break;
    }
  };

  //
  // Initialization
  //
  for (V p = 0; p < num_v; ++p) {
    k1[p] = X[p].size();
    k2[p] = k - k1[p];
    c[p] = 0;
    last_blue[p] = X[p].back();
    if (X[p].size() == (size_t)k) {
      que[0].push({last_blue[p], p});
      type[p] = 0;
      T[k2[p] + 1].insert({last_blue[p], p});
    } else {
      que[0].push({num_v, p});
      que[1].push({k2[p], p});
      type[p] = 1;
      T[k2[p]].insert({last_blue[p], p});
    }
    for (V ri : X[p]) {
      I[ri].insert(p);
    }
  }

  //
  // Main loop
  //
  while ((Xs.size() == k ? *Xs.rbegin() : num_v) > k - 1) {
    // Selection
    V select = -1;
    V argmin = num_v + 1;
    for (int q = 0; q < 2; q++) {
      while (!que[q].empty()) {
        // Remove unnecessary or changed elements from queue
        V top_key = que[q].top().first;
        V top_v = que[q].top().second;

        if (top_key == num_v && Xs.size() < (size_t)k && !centered[top_v]) {
          break;  // Fit to Naive Method
        } else if (removed[top_v] || q != type[top_v]) {
          que[q].pop();
          continue;
        } else if (q == 1 && k2[top_v] != top_key) {
          que[q].pop();
          continue;
        }
        break;
      }
      if (que[q].empty()) continue;
      V v = que[q].top().second;

      set<V> tmp(Xs);
      merge_and_purify(tmp, X[v], k);
      V ec_tmp = tmp.size() == k ? *tmp.rbegin() : num_v;
      if (argmin > ec_tmp || (argmin == ec_tmp && v < select)) {
        argmin = ec_tmp;
        select = v;
      }
    }
    assert(select >= 0);

    //
    // Merge and Remove
    //
    vector<V> delta = merge_and_purify(Xs, X[select], k);
    centers.push_back(select);
    centered[select] = true;
    removed[select] = true;
    remove_multimap_pair(k2[select], {last_blue[select], select});

    //
    // Update about covered elements
    //
    for (V rank_i : delta) {
      covered_rank[rank_i] = true;
      for (V box : I[rank_i]) {
        if (k2[box] + 1 >= k) removed[box] = true;
        if (removed[box]) continue;
        pair<V, V> box_pair = {last_blue[box], box};

        // When the subbox is Type 1 and its last element is covered,
        // it has to be Type2

        if (type[box] == 1) {
          remove_multimap_pair(k2[box], box_pair);
          T[k2[box] + 1].insert(box_pair);
          que[1].push({k2[box] + 1, box});
        } else if (rank_i == last_blue[box]) {
          remove_multimap_pair(k2[box] + 1, box_pair);
          remove_covered_ranks(box);
          if (removed[box]) continue;
          type[box] = 1;
          if (k2[box] + 1 <= k) T[k2[box] + 1].insert({last_blue[box], box});
          que[1].push({k2[box] + 1, box});
        } else {
          remove_multimap_pair(k2[box] + 1, box_pair);
          T[k2[box] + 2].insert(box_pair);
        }

        c[box]++;
        k2[box]++;
      }
    }

    //
    // Update
    //
    int j = 1;
    for (auto it_Xs = Xs.begin(); it_Xs != Xs.end(); ++it_Xs, ++j) {
      if (T[j].empty()) continue;

      vector<pair<V, V>> removing;  // removed elements from T after iteration
      for (auto it = T[j].lower_bound(*it_Xs); it != T[j].end(); ++it) {
        V box = (*it).second;

        auto it_j = it_Xs;
        // always last_blue[box] >= jth_rank
        removing.push_back({last_blue[box], box});
        if (removed[box]) continue;
        if (type[box] == 0) {  // Type1->
          assert(k2[box] + 1 == j);
          k1[box]--;
          k2[box]++;
          if (covered_rank[last_blue[box]]) {
            c[box]--;
            k2[box]--;
            it_j--;
          }
          if (k1[box] == 0 || k2[box] >= k) {
            removed[box] = true;
            continue;
          }
          last_blue[box] = X[box][k1[box] - 1];
        }

        remove_covered_ranks(box);
        if (removed[box]) continue;
        if (last_blue[box] > *it_j) {  // Type1->Type1
          insert_as_type(box, 0);
        } else {  // Type1->Type2
          insert_as_type(box, 1);
        }
      }
      for (auto rm_pair : removing) remove_multimap_pair(j, rm_pair);
    }
  }
}

// Select-Greedily O(n^2*k)
void naive_select_greedily(const G &g, const vector<vector<V>> &X,
                           vector<V> &centers, vector<bool> &centered,
                           const int k) {
  V num_v = g.num_vertices();
  set<V> Xs;
  auto estimated_cardinality = [&](const set<V> &subset) -> double {
    if (subset.size() < (size_t)k) {
      return (k - 1);
    }
    if (subset.size() == (size_t)k) {
      V kth = *subset.rbegin();
      return (double)(k - 1) / kth * g.num_vertices();
    }
    assert(false);
    return (k - 1);
  };

  while (estimated_cardinality(Xs) < max(num_v, k - 1)) {
    V selected_v = -1;
    double argmax = 0.0;
    for (V v = 0; v < num_v; v++) {
      if (centered[v]) continue;
      set<V> tmp(Xs);
      merge_and_purify(tmp, X[v], k);
      double ec_tmp = estimated_cardinality(tmp);
      if (argmax < ec_tmp) {
        argmax = ec_tmp;
        selected_v = v;
      }
    }

    if (selected_v < 0) {
      break;
    }

    centers.push_back(selected_v);
    centered[selected_v] = true;

    merge_and_purify(Xs, X[selected_v], k);
  }
}

vector<V> box_cover_sketch(const G &g, W radius, const int k,
                           const int pass_num, const double aim_coverage) {
  assert(k > 0);

  const V num_v = g.num_vertices();

  vector<V> centers;
  vector<bool> is_covered(num_v, false);
  vector<bool> centered(num_v, false);
  vector<V> rank(num_v);
  vector<V> inv(num_v);
  for (V i = 0; i < num_v; ++i) {
    inv[i] = i;
  }
  for (int pass_trial = 0; pass_trial < pass_num; pass_trial++) {
    //
    // Build-Sketches O((n+m)*rad)
    //
    random_shuffle(inv.begin(), inv.end());
    for (int i = 0; i < num_v; ++i) {
      rank[inv[i]] = i;
    }
    vector<vector<V>> X = build_sketch(g, radius, k, rank, inv, is_covered);

    //
    // Select-Greedily O(n^2*k)
    //
    // naive_select_greedily(g, X, centers, centered, k);
    select_greedily(g, X, centers, centered, k);
    if (coverage(g, centers, radius) >= aim_coverage) break;
  }

  return centers;
}
