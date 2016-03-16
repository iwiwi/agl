#pragma once

DEFINE_int32(gtp_dfs_edge_max, 1000000000, "");
long long gtp_edge_count = 0;
long long gtp_edge_miss = 0;
long long gtp_edge_use = 0;

// ある頂点からdfsをして、貪欲にtree packingを求める
class greedy_treepacking {
  
  class vecE {
  public:
    vecE() {
      idx_ = 0;
    }

    void add(int to) {
      to_.emplace_back(to);
    }

    bool empty() const { return to_.empty(); }
    int size() const { return sz(to_); }
    void advance() {
      idx_++;
      if (sz(to_) == idx_) idx_ = 0;
    }
    const V current() const { return to_[idx_]; }
    void remove_current() { 
      swap(to_[idx_], to_.back());
      to_.pop_back();
      if(sz(to_) == idx_) idx_ = 0;
    }
  public:
    int rem_size_, idx_;
    vector<V> to_;
  };

  void dfs(int v) {
    used_revision_[v] = vertices_revision_;
    auto& to_edges = edges_[v];
    const int rem_edges = min(to_edges.size(), FLAGS_gtp_dfs_edge_max);
    // const int rem_edges = to_edges.size();
    FOR(i, rem_edges) {
      V to = to_edges.current();
      gtp_edge_count++;
      if (used_revision_[to] == vertices_revision_) {
        to_edges.advance();
        gtp_edge_miss++;
      } else {
        to_edges.remove_current();
        inedge_count_[to]++; // in-edgeの本数が増える
        dfs(to);
        gtp_edge_use++;
      }
    }
  }

public:
  greedy_treepacking(const vector<pair<V, V>>& edges, int num_vs) :
    n_(num_vs), edges_(n_), inedge_count_(n_), used_revision_(n_), vertices_revision_(1) {
    for (auto& e : edges) {
      edges_[e.first].add(e.second);
      edges_[e.second].add(e.first);
    }

    for(auto& e : edges_) {
      sort(e.to_.begin(), e.to_.end(), [&edges_ = this->edges_](const V l,const V r) {
        int a = edges_[l].size();
        int b = edges_[r].size();
        return a < b;
      });
    }
  }

  void arborescence_packing(int from) {
    FOR(_, sz(edges_[from])) {
      if(vertices_revision_% 1000 == 0) {
        stringstream ss;
        ss << "revision = " << vertices_revision_ << ", gtp_edge_count = " << gtp_edge_count << ", gtp_edge_use = " << gtp_edge_use;
        JLOG_ADD("greedy_treepacking.progress", ss.str());
      }
      dfs(from);
      vertices_revision_++;
    }
  }

  int inedge_count(int v) const {
    return inedge_count_[v];
  }

private:
  int n_;
  vector<vecE> edges_;
  vector<int> inedge_count_; // dfs tree packingで、その頂点に何本のin-edgeがあるか
  vector<int> used_revision_;
  int vertices_revision_;
};
