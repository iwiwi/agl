#pragma once

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
    void shuffle() {
      random_shuffle(to_.begin(),to_.end(), agl::random);
    }
  private:
    int rem_size_, idx_;
    vector<V> to_;
  };

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

        if(solved_[v] && solved_[to]) {
          edges_[to].add(v); // v,to両方、解が欲しい頂点ではないならば、残余グラフを作る容量で、逆辺を追加する
        }
      }
    }
  }

public:
  greedy_treepacking(const vector<pair<V, V>>& edges, int num_vs, vector<bool>& solved) :
    n_(num_vs), edges_(n_), inedge_count_(n_), used_revision_(n_), vertices_revision_(1), solved_(solved) {
    for (auto& e : edges) {
      edges_[e.first].add(e.second);
      edges_[e.second].add(e.first);
    }
  }

  void arborescence_packing(int from) {
    for(auto& e : edges_) e.shuffle();
    FOR(_, sz(edges_[from])) {
      dfs(from);
      vertices_revision_++;
    }
  }

  int inedge_count(int v) const {
    return inedge_count_[v];
  }

  void set_solved(int v) {
    solved_[v] = true;
  }

private:
  const int n_;
  vector<vecE> edges_;
  vector<int> inedge_count_; // dfs tree packingで、その頂点に何本のin-edgeがあるか
  vector<int> used_revision_;
  int vertices_revision_;
  vector<bool>& solved_;
};
