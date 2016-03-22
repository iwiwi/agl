#pragma once
#include "graph/graph.h"
#include "graph/graph_index_interface.h"
#include <set>
#include <queue>

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

  static const W W_INF = 100;
  struct index_t {
    W bpspt_d[kNumBitParallelRoots];
    uint64_t bpspt_s[kNumBitParallelRoots][2];
    std::map<V, W> spm;
    void update(V v, W d) {
      if (spm.count(v) && spm[v] <= d) return;
      spm[v] = d;
    }
  };

  std::vector<index_t> idx[2];
  std::vector<std::vector<V>> adj[2];
  std::vector<V> rank;
  std::vector<V> inv;

 private:
  void pruned_bfs(V root, int direction);
  void resume_pbfs(V v_from, V v_to, W d_ft, int direction);
  W query_distance(V v_from, V v_to, int direction);
};

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::construct(
    const G &g) {
  V num_v = g.num_vertices();
  assert(num_v >= 3);
  adj[0].resize(num_v), adj[1].resize(num_v);
  idx[0].resize(num_v), idx[1].resize(num_v);

  for (const auto &p : g.edge_list()) {
    adj[0][p.first].push_back(p.second);
    adj[1][p.second].push_back(p.first);
  }
  for (V v = 0; v < num_v; ++v) {
    std::sort(adj[0][v].begin(), adj[0][v].end());
    std::sort(adj[1][v].begin(), adj[1][v].end());
  }

  rank.resize(num_v);
  inv.resize(num_v);
  {
    std::vector<std::pair<double, V>> decreasing_order(num_v);
    for (V i = 0; i < num_v; ++i) {
      double deg = adj[0][i].size() + adj[1][i].size();
      double t = (double)agl::random(num_v) / num_v;
      decreasing_order[i] = {deg + t, i};
    }
    std::sort(decreasing_order.rbegin(), decreasing_order.rend());
    for (int i = 0; i < num_v; ++i) inv[i] = decreasing_order[i].second;
    for (int i = 0; i < num_v; ++i) rank[inv[i]] = i;
  }

  std::vector<bool> used(num_v, false);
  {
    V root = 0;
    for (int bp_i = 0; bp_i < kNumBitParallelRoots; ++bp_i) {
      // Select Root
      while (root < num_v && used[root]) root++;
      if (root == num_v) {
        for (V v = 0; v < num_v; ++v) {
          idx[0][v].bpspt_d[bp_i] = W_INF;
          idx[1][v].bpspt_d[bp_i] = W_INF;
        }
        continue;
      }
      used[root] = true;

      // Select Roots
      std::set<V> neighbor_ranks;
      for (V v : adj[0][root]) neighbor_ranks.insert(rank[v]);
      for (V v : adj[1][root]) neighbor_ranks.insert(rank[v]);

      std::vector<V> selected_roots;
      for (V r : neighbor_ranks) {
        V v = inv[r];
        if (used[v]) continue;
        used[v] = true;
        selected_roots.push_back(v);
        if (selected_roots.size() == 64) break;
      }

      for (int direction = 0; direction < 2; ++direction) {
        int another = direction ^ 1;

        std::vector<W> tmp_d(num_v, W_INF);
        std::vector<std::pair<uint64_t, uint64_t>> tmp_s(num_v, {0, 0});
        std::queue<std::pair<V, W>> que;

        que.emplace(root, 0);
        tmp_d[root] = 0;

        for (size_t i = 0; i < selected_roots.size(); ++i) {
          V v = selected_roots[i];
          if (!std::binary_search(adj[direction][root].begin(),
                                  adj[direction][root].end(), v))
            continue;
          que.emplace(v, 1);
          tmp_d[v] = 1;
          tmp_s[v].first = 1ULL << i;
        }

        for (W d = 0; !que.empty(); ++d) {
          std::vector<std::pair<V, V>> sibling_es;
          std::vector<std::pair<V, V>> child_es;
          while (!que.empty() && que.front().second == d) {
            V v = que.front().first;
            que.pop();

            for (V tv : adj[direction][v]) {
              W td = d + 1;
              if (d == tmp_d[tv]) {
                if (v >= tv) continue;
                sibling_es.emplace_back(v, tv);
              } else if (d < tmp_d[tv]) {
                if (tmp_d[tv] == W_INF) {
                  que.emplace(tv, td);
                  tmp_d[tv] = td;
                }
                child_es.emplace_back(v, tv);
              }
            }
          }

          for (auto &p : sibling_es) {
            V v, w;
            std::tie(v, w) = p;
            tmp_s[v].second |= tmp_s[w].first;
            tmp_s[w].second |= tmp_s[v].first;
          }
          for (auto &p : child_es) {
            V v, c;
            std::tie(v, c) = p;
            tmp_s[c].first |= tmp_s[v].first;
            tmp_s[c].second |= tmp_s[v].second;
          }
        }

        for (V v = 0; v < num_v; ++v) {
          idx[another][v].bpspt_d[bp_i] = tmp_d[v];
          idx[another][v].bpspt_s[bp_i][0] = tmp_s[v].first;
          idx[another][v].bpspt_s[bp_i][1] = tmp_s[v].second & ~tmp_s[v].first;
        }
      }
    }
  }

  // Pruned labelling
  for (const auto &v : inv) {
    if (used[v]) continue;
    pruned_bfs(v, 0);
    pruned_bfs(v, 1);
  }
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::pruned_bfs(
    V root, int direction) {
  V num_v = adj[0].size();
  std::queue<V> que;
  que.push(root);
  std::vector<W> P(num_v, W_INF);
  P[root] = 0;

  std::vector<std::pair<V, W>> tmp_idx;

  while (!que.empty()) {
    V u = que.front();
    que.pop();
    if (u != root && query_distance(root, u, direction) <= P[u]) continue;
    tmp_idx.emplace_back(u, P[u]);

    for (const auto &w : adj[direction][u]) {
      if (P[w] < W_INF) continue;
      P[w] = P[u] + 1;
      que.push(w);
    }
  }

  int another = direction ^ 1;
  for (const auto &p : tmp_idx) idx[another][p.first].update(root, p.second);
}

template <size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance(
    const G &g, V v_from, V v_to) {
  return query_distance(v_from, v_to, 0);
}

template <size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance(
    V v_from, V v_to, int direction) {
  int another = direction ^ 1;
  assert(v_from < adj[0].size() && v_to < adj[0].size());
  if (v_from == v_to) return 0;

  W d = W_INF;
  const index_t &idx_from = idx[direction][v_from];
  const index_t &idx_to = idx[another][v_to];

  // TODO bit-parallel
  for (int bp_i = 0; bp_i < kNumBitParallelRoots; ++bp_i) {
    W td = idx_from.bpspt_d[bp_i] + idx_to.bpspt_d[bp_i];
    if (td - 2 > d) continue;
    if (idx_from.bpspt_s[bp_i][0] & idx_to.bpspt_s[bp_i][0]) {
      td -= 2;
    } else if ((idx_from.bpspt_s[bp_i][1] & idx_to.bpspt_s[bp_i][0]) |
               (idx_from.bpspt_s[bp_i][0] & idx_to.bpspt_s[bp_i][1])) {
      td -= 1;
    }
    if (td < d) d = td;
  }

  for (auto &p : idx_from.spm)
    if (p.first == v_to && d > p.second) d = p.second;
  for (auto &p : idx_to.spm)
    if (p.first == v_from && d > p.second) d = p.second;

  for (auto i1 = idx_from.spm.begin(), i2 = idx_to.spm.begin();
       i1 != idx_from.spm.end() && i2 != idx_to.spm.end();) {
    V v1 = i1->first;
    V v2 = i2->first;
    if (v1 == v2) {
      W td = i1->second + i2->second;
      if (td < d) d = td;
      i1++;
      i2++;
    } else {
      if (v1 < v2) i1++;
      if (v1 > v2) i2++;
    }
  }
  if (d >= W_INF - 2) d = W_INF;
  return d;
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::resume_pbfs(
    V neighbor_a, V v_b, W d_nab, int direction) {
  int another = direction ^ 1;

  std::queue<std::pair<V, W>> que;
  que.emplace(v_b, d_nab);
  while (!que.empty()) {
    V v;
    W d;
    std::tie(v, d) = que.front();
    que.pop();
    if (query_distance(neighbor_a, v, direction) <= d) continue;
    idx[another][v].update(neighbor_a, d);
    for (const auto &w : adj[direction][v]) que.emplace(w, d + 1);
  }
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::add_edge(
    const G &g, V v_a, const E &e) {
  V v_b = to(e);
  assert(v_a < g.num_vertices() && v_b < g.num_vertices());
  if (query_distance(g, v_a, v_b) <= 1) return;
  adj[0][v_a].push_back(v_b);
  adj[1][v_b].push_back(v_a);
  std::sort(adj[0][v_a].begin(), adj[0][v_a].end());
  std::sort(adj[1][v_b].begin(), adj[1][v_b].end());

  std::vector<std::tuple<V, W, int>> tmp;
  for (const auto &p : idx[1][v_a].spm) {
    V r = rank[p.first];
    V d = p.second;
    tmp.emplace_back(r, d, 0);
  }
  for (const auto &p : idx[0][v_b].spm) {
    V r = rank[p.first];
    V d = p.second;
    tmp.emplace_back(r, d, 1);
  }

  std::sort(tmp.begin(), tmp.end());

  for (const auto &q : tmp) {
    V r, t;
    W d;
    std::tie(r, d, t) = q;
    if (t == 0) {
      V neighbor_a = inv[r];
      resume_pbfs(neighbor_a, v_b, d + 1, 0);
    }
    if (t == 1) {
      V neighbor_b = inv[r];
      resume_pbfs(neighbor_b, v_a, d + 1, 1);
    }
  }
}
}  // namespace agl
