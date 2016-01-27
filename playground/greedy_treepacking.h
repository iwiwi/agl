#pragma once

// ある頂点からdfsをして、貪欲にtree packingを求める
class greedy_treepacking {
  
  class vecE {
  public:
    void add(int to) {
      to_.emplace_back(to);
    }

    void reset() {
      rem_size_ = sz(to_);
      idx_ = 0;
    }

    bool empty() const { return rem_size_ == 0; }
    int size() const { return rem_size_; }
    void advance() {
      idx_++;
      if (rem_size_ == idx_) idx_ = 0;
    }
    const V current() const { return to_[idx_]; }
    void remove_current() { 
      swap(to_[rem_size_ - 1], to_[idx_]);
      rem_size_--;
      if (idx_ == rem_size_) idx_ = 0;
    }

  private:
    int rem_size_, idx_;
    vector<V> to_;
  };

  void reset_graph() {
    for (auto& ev : edges_) ev.reset();
    fill(inedge_count_.begin(), inedge_count_.end(), 0);

    if (vertices_revision_ >= numeric_limits<decltype(vertices_revision_)>::max() / 2) {
      vertices_revision_ = 1;
      fill(used_revision_.begin(), used_revision_.end(), 0);
    }
  }

  void dfs(int v) {
    used_revision_[v] = vertices_revision_;
    auto& to_edges = edges_[v];
    const int rem_edges = to_edges.size();
    FOR(i, rem_edges) {
      V to = to_edges.current();
      if (used_revision_[to] == vertices_revision_) {
        to_edges.advance();
        continue;
      } else {
        to_edges.remove_current();
        inedge_count_[to]++; // in-edgeの本数が増える
        dfs(to);
      }
    }
  }

public:
  greedy_treepacking(vector<pair<V, V>>& edges, int num_vs) :
    n_(num_vs), edges_(n_), inedge_count_(n_), used_revision_(n_), vertices_revision_(1) {
    for (auto& e : edges) {
      edges_[e.first].add(e.second);
      edges_[e.second].add(e.first);
    }
  }

  void arborescence_packing(int from) {
    reset_graph();
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
  vector<vecE> edges_;
  vector<int> inedge_count_; // dfs tree packingで、その頂点に何本のin-edgeがあるか
  vector<int> used_revision_;
  int vertices_revision_;
};
