#pragma once
#include "graph/graph.h"
#include "graph/graph_index_interface.h"

namespace agl {
template <size_t kNumBitParallelRoots = 16>
class dynamic_pruned_landmark_labeling
    : public dynamic_graph_index_interface<G>,
      public distance_query_interface<G> {
 public:
  virtual ~dynamic_pruned_landmark_labeling() {}
  virtual void construct(const G &g) override;
  virtual W query_distance(const G &g, V v_from, V v_to) override;
  virtual void add_edge(const G &g, V v_from, const E &e) override;
  virtual void remove_edge(const G &g, V v_from, V v_to) override {
    assert(false);
  }
  virtual void add_vertices(const G &g, V old_num_vertices) override {
    assert(false);
  }
  virtual void remove_vertices(const G &, V old_num_vertices) override {
    assert(false);
  }

  const uint8_t D_INF = 100;
  struct index_t {
    uint8_t bpspt_d[kNumBitParallelRoots];
    uint64_t bpspt_s[kNumBitParallelRoots][2];
    std::vector<V> spt_v;
    std::vector<uint8_t> spt_d;
    void update(V v, uint8_t d) {
      if (spt_v.empty() || spt_v.back() < v) {
        spt_v.emplace_back(v);
        spt_d.emplace_back(d);
        return;
      }
      size_t i;
      for (i = 0; i < spt_v.size(); ++i) {
        if (spt_v[i] >= v) break;
      }
      assert(i < spt_v.size());
      if (spt_v[i] == v) {
        if (spt_d[i] > d) spt_d[i] = d;
        return;
      }
      assert(spt_v[i] > v);
      spt_v.push_back(spt_v.back());
      spt_d.push_back(spt_d.back());
      for (size_t j = spt_v.size() - 1; j > i; --j) {
        spt_v[j] = spt_v[j - 1];
        spt_d[j] = spt_d[j - 1];
      }
      spt_v[i] = v;
      spt_d[i] = d;
    }
    size_t size() const { return spt_v.size(); }
  };

  V num_v;
  std::vector<index_t> idx[2];
  std::vector<std::vector<V>> adj[2];
  std::vector<V> rank;
  std::vector<V> inv;

  size_t total_label_num() {
    size_t sum = 0;
    for (int i = 0; i < 2; ++i)
      for (const index_t &j : idx[i]) sum += j.size();
    return sum;
  }

  // private:
  void load_graph(const G &g);
  void pruned_bfs(V root, int direction, const std::vector<bool> &used);
  void partial_bfs(V v_from, V v_to, uint8_t d_ft, int direction);
  uint8_t distance_less(V v_from, V v_to, int direction, uint8_t upper_limit);
  void bit_parallel_bfs(std::vector<bool> &used, const size_t num_e);
  void partial_bp_bfs(int bp_i, V v_from, int direction);
  // Reusable containers
  std::vector<uint8_t> bfs_dist;
  std::vector<V> bfs_que;
  V bp_roots[64];
};

/**
* Construct the data structure from the given graph.
* \param g graph
*/
template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::construct(
    const G &g) {
  num_v = g.num_vertices();
  size_t num_e = g.num_edges();
  assert(num_v >= 3);

  // Initialize
  double timer = -get_current_time_sec();
  load_graph(g);
  timer += get_current_time_sec();
  std::cerr << timer << " sec load" << std::endl;

  std::vector<bool> used(num_v, false);

  // Bit-Parallel Labeling
  timer = -get_current_time_sec();
  bit_parallel_bfs(used, num_e);
  timer += get_current_time_sec();
  std::cerr << timer << " sec bit-parallel" << std::endl;

  timer = -get_current_time_sec();
  // Pruned labelling
  for (V root = 0; root < num_v; ++root) {
    if (used[root]) continue;
    pruned_bfs(root, 0, used);
    pruned_bfs(root, 1, used);
    used[root] = true;
  }
  timer += get_current_time_sec();
  std::cerr << timer << " sec pruned-bfs" << std::endl;
}

/**
* Construct the Bit-Parallel label
* \param used vector to record the vertices used as the root of BFS
* \param num_e the number of edges of the given graph
*/
template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::bit_parallel_bfs(
    std::vector<bool> &used, const size_t num_e) {
  V root = 0;
  std::vector<std::pair<uint64_t, uint64_t>> tmp_s(num_v, {0, 0});
  std::vector<std::pair<V, V>> sibling_es(num_e);
  std::vector<std::pair<V, V>> child_es(num_e);
  for (int bp_i = 0; bp_i < kNumBitParallelRoots; ++bp_i) {
    // Select Root
    while (root < num_v && used[root]) root++;
    if (root == num_v) {
      for (V v = 0; v < num_v; ++v) {
        idx[0][v].bpspt_d[bp_i] = D_INF;
        idx[1][v].bpspt_d[bp_i] = D_INF;
      }
      continue;
    }
    used[root] = true;

    // Select Roots
    V selected_num = 0;
    for (int o = 0, i = 0;
         o < adj[0][root].size() && i < adj[1][root].size();) {
      V vo = adj[0][root][o], vi = adj[1][root][i];
      if (vo == vi) {
        if (!used[vo]) {
          used[vo] = true;
          bp_roots[selected_num++] = vo;
          if (selected_num == 64) break;
        }
        o++;
        i++;
      } else {
        if (vo > vi) i++;
        if (vo < vi) o++;
      }
    }

    for (int direction = 0; direction < 2; ++direction) {
      int another = direction ^ 1;

      int q_head = 0, q_tail = 0;

      bfs_que[q_tail++] = root;
      bfs_dist[root] = 0;

      for (size_t i = 0; i < selected_num; ++i) {
        V v = bp_roots[i];
        bfs_que[q_tail++] = v;
        bfs_dist[v] = 1;
        tmp_s[v].first = 1ULL << i;
      }

      for (uint8_t d = 0; q_head < q_tail; ++d) {
        size_t num_sibling_es = 0, num_child_es = 0;
        while (q_head < q_tail && bfs_dist[bfs_que[q_head]] == d) {
          V v = bfs_que[q_head++];
          for (const V &tv : adj[direction][v]) {
            uint8_t td = d + 1;
            if (d > bfs_dist[tv]) continue;
            if (d == bfs_dist[tv] && v < tv) {
              sibling_es[num_sibling_es].first = v;
              sibling_es[num_sibling_es].second = tv;
              num_sibling_es++;
            }
            if (d < bfs_dist[tv]) {
              child_es[num_child_es].first = v;
              child_es[num_child_es].second = tv;
              num_child_es++;
            }
            if (bfs_dist[tv] == D_INF) {
              bfs_dist[tv] = td;
              bfs_que[q_tail++] = tv;
            }
          }
        }

        for (int p = 0; p < num_sibling_es; ++p) {
          V v, w;
          std::tie(v, w) = sibling_es[p];
          tmp_s[v].second |= tmp_s[w].first;
          tmp_s[w].second |= tmp_s[v].first;
        }
        for (int p = 0; p < num_child_es; ++p) {
          V v, c;
          std::tie(v, c) = child_es[p];
          tmp_s[c].first |= tmp_s[v].first;
          tmp_s[c].second |= tmp_s[v].second;
        }
      }

      for (V v = 0; v < num_v; ++v) {
        idx[another][v].bpspt_d[bp_i] = bfs_dist[v];
        idx[another][v].bpspt_s[bp_i][0] = tmp_s[v].first;
        idx[another][v].bpspt_s[bp_i][1] = tmp_s[v].second & ~tmp_s[v].first;
        tmp_s[v].first = 0, tmp_s[v].second = 0;
      }
      for (int i = 0; i < q_tail; ++i) bfs_dist[bfs_que[i]] = D_INF;
    }
  }
}

/**
* Load and relabel the given graph
* \param g graph
*/
template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::load_graph(
    const G &g) {
  rank.resize(num_v);
  inv.resize(num_v);
  {
    std::vector<V> deg(num_v, 0);
    for (const auto &p : g.edge_list())
      if (p.first != p.second) deg[p.first]++, deg[p.second]++;

    std::vector<std::pair<double, V>> sorting_v(num_v);
    for (V v = 0; v < num_v; ++v) {
      double t = (double)agl::random(num_v) / num_v;
      sorting_v[v] = std::make_pair(deg[v] + t, v);
    }
    std::sort(sorting_v.rbegin(), sorting_v.rend());
    for (int i = 0; i < num_v; ++i) inv[i] = sorting_v[i].second;
    for (int i = 0; i < num_v; ++i) rank[inv[i]] = i;
  }

  idx[0].resize(num_v), idx[1].resize(num_v);
  adj[0].resize(num_v), adj[1].resize(num_v);
  for (const auto &p : g.edge_list()) {
    V v = p.first, u = p.second;
    if (v == u) continue;
    adj[0][rank[v]].push_back(rank[u]);
    adj[1][rank[u]].push_back(rank[v]);
  }
  for (V v = 0; v < num_v; ++v) {
    std::sort(adj[0][v].begin(), adj[0][v].end());
    std::sort(adj[1][v].begin(), adj[1][v].end());
  }
  bfs_dist.assign(num_v, D_INF);
  bfs_que.resize(num_v);
}

/**
* Construct 2-Hop label by pruned BFS
* \param root the root vertex of pruned BFS
* \param direction the direction to do BFS
* \param used vector to record the vertices used as the root of BFS
*/
template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::pruned_bfs(
    V root, int direction, const std::vector<bool> &used) {
  int another = direction ^ 1;
  int q_head = 0, q_tail = 0;
  bfs_que[q_tail++] = root;
  bfs_dist[root] = 0;

  while (q_head < q_tail) {
    V u = bfs_que[q_head++];
    if (u != root &&
        distance_less(root, u, direction, bfs_dist[u]) <= bfs_dist[u])
      continue;
    idx[another][u].update(root, bfs_dist[u]);

    for (const auto &w : adj[direction][u]) {
      if (used[w]) continue;
      if (bfs_dist[w] < D_INF) continue;
      bfs_dist[w] = bfs_dist[u] + 1;
      bfs_que[q_tail++] = w;
    }
  }
  for (int i = 0; i < q_tail; ++i) {
    bfs_dist[bfs_que[i]] = D_INF;
  }
}

/**
* Return the distance of shortest path between two vertices.
* \param g graph
* \param v_from the starting vertex of shortest path
* \param v_to the end vertex of shortest path
*/
template <size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance(
    const G &g, V v_from, V v_to) {
  assert(v_from >= 0 && v_to >= 0);
  assert(v_from < rank.size() && v_to < rank.size());
  if (v_from == v_to) return 0;
  v_from = rank[v_from], v_to = rank[v_to];
  return distance_less(v_from, v_to, 0, 0);
}

/**
* Return the greater either the distance of shortest path between two vertices
* or upper_limit
* \param v_from the starting vertex of shortest path
* \param v_to the end vertex of shortest path
* \param direction the direction of shortest path
* \param upper_limit When it is found out that the distance of the shortest path
* is smaller than upper_limit, this function returns the calculationg distance
* which is equal to or shorter than upper_limit.
*/
template <size_t kNumBitParallelRoots>
uint8_t dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::distance_less(
    V v_from, V v_to, int direction, uint8_t upper_limit) {
  int another = direction ^ 1;

  uint8_t d = D_INF;
  const index_t &idx_from = idx[direction][v_from];
  const index_t &idx_to = idx[another][v_to];

  for (int bp_i = 0; bp_i < kNumBitParallelRoots; ++bp_i) {
    uint8_t td = idx_from.bpspt_d[bp_i] + idx_to.bpspt_d[bp_i];
    if (td - 2 >= d) continue;
    if (idx_from.bpspt_s[bp_i][0] & idx_to.bpspt_s[bp_i][0]) {
      td -= 2;
    } else if ((idx_from.bpspt_s[bp_i][1] & idx_to.bpspt_s[bp_i][0]) |
               (idx_from.bpspt_s[bp_i][0] & idx_to.bpspt_s[bp_i][1])) {
      td -= 1;
    }
    if (td < d) {
      d = td;
      if (d <= upper_limit) return d;
    }
  }

  for (size_t i1 = 0, i2 = 0; i1 < idx_from.size() && i2 < idx_to.size();) {
    V v1 = idx_from.spt_v[i1];
    V v2 = idx_to.spt_v[i2];
    if (v1 == v2) {
      uint8_t td = idx_from.spt_d[i1] + idx_to.spt_d[i2];
      if (td < d) {
        d = td;
        if (d <= upper_limit) return d;
      }
      i1++;
      i2++;
    } else {
      if (v1 == v_to && idx_from.spt_d[i1] <= d) return idx_from.spt_d[i1];
      if (v2 == v_from && idx_to.spt_d[i2] <= d) return idx_to.spt_d[i2];
      if (v1 < v2) i1++;
      if (v1 > v2) i2++;
    }
  }
  if (d >= D_INF - 2) d = D_INF;
  return d;
}

/**
* Update 2-Hop label when the new edge is added.
* \param v_from the root vertex of partial BFS
* \param v_end one of the end vertex of new edge
* \param base_d the distance between v_from to v_end
* \param direction direction
*/
template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::partial_bfs(
    V v_from, V v_end, uint8_t base_d, int direction) {
  int another = direction ^ 1;

  int q_head = 0, q_tail = 0;
  bfs_que[q_tail++] = v_end;
  bfs_dist[v_end] = base_d;
  while (q_head < q_tail) {
    V v = bfs_que[q_head++];
    uint8_t d = bfs_dist[v];
    if (distance_less(v_from, v, direction, d) <= d) continue;
    idx[another][v].update(v_from, d);
    for (const auto &w : adj[direction][v]) {
      if (bfs_dist[w] < D_INF) continue;
      bfs_dist[w] = d + 1;
      bfs_que[q_tail++] = w;
    }
  }
  for (int i = 0; i < q_tail; ++i) {
    bfs_dist[bfs_que[i]] = D_INF;
  }
}

/**
* Update bit-parallel label when the new edge is added.
* \param bp_i index of bit-parallel
* \param v_from the root vertex of partial BFS
* \param direction direction
*/
template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::partial_bp_bfs(
    int bp_i, V v_from, int direction) {
  V another = direction ^ 1;
  const index_t &idx_from = idx[another][v_from];
  const uint8_t base_d = idx_from.bpspt_d[bp_i];
  if (base_d == D_INF) return;

  int q_head = 0, q_tail = 0;
  bfs_que[q_tail++] = v_from;
  bfs_dist[v_from] = base_d;
  for (uint8_t d = base_d; q_head < q_tail; ++d) {
    int old_head = q_head;
    while (q_head < q_tail && bfs_dist[bfs_que[q_head]] == d) {
      V v = bfs_que[q_head++];
      const index_t &idx_v = idx[another][v];

      for (V tv : adj[direction][v]) {
        index_t &idx_tv = idx[another][tv];
        if (d == idx_tv.bpspt_d[bp_i])
          idx_tv.bpspt_s[bp_i][1] |= idx_v.bpspt_s[bp_i][0];
        if (d + 1 < idx_tv.bpspt_d[bp_i]) {
          idx_tv.bpspt_s[bp_i][1] =
              (idx_tv.bpspt_d[bp_i] == d + 2 ? idx_tv.bpspt_s[bp_i][0] : 0);

          idx_tv.bpspt_s[bp_i][0] = 0;
          idx_tv.bpspt_d[bp_i] = d + 1;
          bfs_dist[tv] = d + 1;
          bfs_que[q_tail++] = tv;
        }
      }
    }

    for (int i = old_head; i < q_head; ++i) {
      V v = bfs_que[i];
      const index_t &idx_v = idx[another][v];
      for (V tv : adj[direction][v]) {
        index_t &idx_tv = idx[another][tv];
        if (idx_tv.bpspt_d[bp_i] == d + 1) {
          // Set propagation (2)
          idx_tv.bpspt_s[bp_i][0] |= idx_v.bpspt_s[bp_i][0];
          idx_tv.bpspt_s[bp_i][1] |= idx_v.bpspt_s[bp_i][1];
        }
      }
    }
  }
  for (int i = 0; i < q_tail; ++i) {
    bfs_dist[bfs_que[i]] = D_INF;
  }
}

/**
* Add new edge and update the labels
* \param g graph
* \param v_a the vertex from which new edge starts
* \param e new edge
*/
template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::add_edge(
    const G &g, V v_a, const E &e) {
  V v_b = to(e);
  assert(v_a >= 0 && v_b >= 0);
  assert(v_a < num_v && v_b < num_v);
  v_a = rank[v_a], v_b = rank[v_b];
  if (std::binary_search(adj[0][v_a].begin(), adj[0][v_a].end(), v_b)) return;

  adj[0][v_a].push_back(v_b);
  adj[1][v_b].push_back(v_a);
  std::sort(adj[0][v_a].begin(), adj[0][v_a].end());
  std::sort(adj[1][v_b].begin(), adj[1][v_b].end());

  for (int bp_i = 0; bp_i < kNumBitParallelRoots; ++bp_i) {
    partial_bp_bfs(bp_i, v_a, 0);
    partial_bp_bfs(bp_i, v_b, 1);
  }

  index_t &idx_a = idx[1][v_a];
  index_t &idx_b = idx[0][v_b];
  for (int ia = 0, ib = 0; ia < idx_a.size() || ib < idx_b.size();) {
    V r_a = ia < idx_a.size() ? idx_a.spt_v[ia] : num_v;
    V r_b = ib < idx_b.size() ? idx_b.spt_v[ib] : num_v;
    if (r_a == r_b) {
      partial_bfs(r_a, v_b, idx_a.spt_d[ia] + 1, 0);
      partial_bfs(r_b, v_a, idx_b.spt_d[ib] + 1, 1);
      ia++;
      ib++;
    } else if (r_a > r_b) {
      partial_bfs(r_b, v_a, idx_b.spt_d[ib] + 1, 1);
      ib++;
    } else {
      partial_bfs(r_a, v_b, idx_a.spt_d[ia] + 1, 0);
      ia++;
    }
  }
}
}  // namespace agl
