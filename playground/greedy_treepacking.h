#pragma once

// ある頂点からdfsをして、貪欲にtree packingを求める
class greedy_treepacking {
  class E {
  public:
    E(int to) : to(to) {
      reset();
    }

    void reset() {
      used_revision = 0;
    }

    bool enable(int revision) const {
      return used_revision != revision;
    }

    void mark_used(int revision) {
      used_revision = revision;
    }

    const int to;
  private:
    int used_revision;
  };

  void reset_graph() {
    for (auto& to_edges : edges_) {
      for (auto& e : to_edges) e.reset();
    }
    fill(used_revision_.begin(), used_revision_.end(), 0);
    vertices_revision_ = 1;
    edges_revision_ = 1;
  }

  void reset_reivision() {
    vertices_revision_++;
    edges_revision_++;
    if (vertices_revision_ >= numeric_limits<decltype(vertices_revision_)>::max() / 2) {
      reset_graph();
    }
  }

  void increase_inedge_count(int v) {
    inedge_count_[v]++;
  }

  void dfs(int v) {
    used_revision_[v] = vertices_revision_;
    auto& to_edges = edges_[v];
    FOR(i, sz(to_edges)) {
      auto& e = to_edges[dfs_iter_[v]++];
      if (dfs_iter_[v] == sz(to_edges)) dfs_iter_[v] = 0;

      if (!e.enable(edges_revision_)) continue;
      if (used_revision_[e.to] == vertices_revision_) continue;

      e.mark_used(edges_revision_);
      increase_inedge_count(e.to); // in-edgeの本数が増える
      dfs(e.to);
    }
  }

public:
  greedy_treepacking(vector<pair<V, V>>& edges, int num_vs) :
    n_(num_vs), edges_(n_), dfs_iter_(n_), used_revision_(n_), inedge_count_(n_),
    vertices_revision_(1), edges_revision_(1) {
    for (auto& e : edges) {
      edges_[e.first].emplace_back(e.second);
      edges_[e.second].emplace_back(e.first);
    }
    // reset_graph();
  }

  void arborescence_packing(int from) {
    reset_reivision();
    fill(inedge_count_.begin(), inedge_count_.end(), 0);
    FOR(_, sz(edges_[from])) {
      dfs(from);
      vertices_revision_++;
    }
  }

  int inedge_count(int v) const {
    return inedge_count_[v];
  }

private:
  const int n_;
  vector<vector<E>> edges_;

  //ある頂点から、どの辺まで使ったかを持つ
  //dinicとは違いオーダーは落ちないがhuristicsとしては効くはず
  vector<int> dfs_iter_;
  vector<int> used_revision_; // dfsでその頂点を訪れたか
  vector<int> inedge_count_; // dfs tree packingで、その頂点に何本のin-edgeがあるか
  int vertices_revision_, edges_revision_; //頂点, 辺のrevision
};
