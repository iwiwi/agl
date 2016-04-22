#include "special_purpose_reachability_index.h"
#include "common.h"
#include <cassert>
#include <queue>
using namespace std;

namespace betweenness_centrality {
  namespace special_purpose_reachability_index {
    inline void SetBit(int &number, int pos, bool bit){
      number &= ~(1 << pos);
      number |= (int(bit) << pos);
    }

    inline bool GetBit(int number, int pos){
      return (number >> pos) & 1;
    }
  
    DynamicSPT::DynamicSPT(int r, vector<vector<int> >  *fadj, vector<vector<int> >  *badj) : root(r), fadj(fadj), badj(badj)
    {
      CHECK(ValidNode(r));
      curr_dist.resize(fadj->size(), INF);
      next_dist.resize(fadj->size(), -1);
      Build();
    }

    void DynamicSPT::ChangeRoot(int new_root){
      vector<int> old_dist(curr_dist);
      fill(curr_dist.begin(), curr_dist.end(), INF);
      fill(next_dist.begin(), next_dist.end(), -1);
      root = new_root;
      CHECK(ValidNode(new_root));
      Build();
    }
  
    void DynamicSPT::InsertEdge(int u, int v){
      chg_nodes.clear();
    
      queue<int> que;
      CHECK(ValidNode(u) && ValidNode(v));
    
      if (curr_dist[v] > curr_dist[u] + 1){
        curr_dist[v] = curr_dist[u] + 1;
        que.push(v);
        chg_nodes.push_back(v);
      }
    
      while (!que.empty()){
        int v = que.front(); que.pop();
        for (int w : fadj->at(v)){
          if (curr_dist[w] > curr_dist[v] + 1){
            curr_dist[w] = curr_dist[v] + 1;
            que.push(w);
            chg_nodes.push_back(w);
          }
        }
      }
      // cerr << "HOGE" << endl;
      // assert(chg_nodes.empty());
    }

    void DynamicSPT::DeleteEdge(int u, int v){
      // Frigioni 2000
      // colorの代わりにnext_distの値を見て色々判定している
      chg_nodes.clear();
      CHECK(ValidNode(u) && ValidNode(v));
    
      if (curr_dist[v] == curr_dist[u] + 1){
        vector<int> upd_nodes;
        CollectChanges({v}, upd_nodes);
        FixChanges(upd_nodes);
        upd_nodes.clear();
      }
    }

    void DynamicSPT::InsertNode(int u){
      CHECK(ValidNode(u));
      while (curr_dist.size() < fadj->size()) curr_dist.push_back(INF);
      while (next_dist.size() < fadj->size()) next_dist.push_back(-1);
    }

    void DynamicSPT::DeleteNode(int u, const vector<int> &u_out, const vector<int> &){
      CHECK(ValidNode(u) && fadj->at(u).empty() && badj->at(u).empty());
    
      if (curr_dist[u] < INF){
        vector<int> start_nodes;
        vector<int> upd_nodes;
        for (int v : u_out){
          if (curr_dist[u] + 1 == curr_dist[v]) start_nodes.push_back(v);
        }
        curr_dist[u] = INF;      
        chg_nodes.push_back(u);
        CollectChanges(start_nodes, upd_nodes);
        FixChanges(upd_nodes);
        CHECK(curr_dist[u] == INF);
      }
    }
    
    void DynamicSPT::Build(){
      queue<int> que;
      curr_dist[root] = 0;
      que.push(root);
      while (!que.empty()){
        int v = que.front(); que.pop();
        for (int w : fadj->at(v)){
          if (curr_dist[w] == INF){
            curr_dist[w] = curr_dist[v] + 1;
            que.push(w);
          }
        }
      }
    }
  
    int DynamicSPT::FindParent(int v){
      // Find the neighbor s of v with minimal distance to s
      for (int w : badj->at(v)){
        // only see neighbors whose distance are not changed.
        bool not_red = next_dist[w] == -1 || next_dist[w] == curr_dist[w];
        if (not_red && curr_dist[w] + 1 == curr_dist[v]){
          return w;
        }
      }
      return -1;
    }

    void DynamicSPT::CollectChanges(const vector<int> &start_nodes, vector<int> &upd_nodes){
      queue<int> que;
      for (int v : start_nodes){
        if (FindParent(v) == -1){
          que.push(v);
          next_dist[v] = INF;
          upd_nodes.push_back(v);
        }
      }
      
      while (!que.empty()){
        int v  = que.front(); que.pop(); 
        for (int w : fadj->at(v)){
          if (next_dist[w] != -1 || curr_dist[w] != curr_dist[v] + 1) continue;

          next_dist[w] = curr_dist[w];
          upd_nodes.push_back(w);
          if (FindParent(w) == -1){
            next_dist[w] = INF;
            que.push(w);
          }
        }
      }
    }

    void DynamicSPT::FixChanges(const vector<int> &upd_nodes){
      // step 3.a
      typedef pair<int, int> P;
      priority_queue<P, vector<P>, greater<P> > que;
    
      for (int v : upd_nodes){
        if (next_dist[v] > curr_dist[v]){
          for (int w : badj->at(v)){
            if (next_dist[w] == -1 || next_dist[w] == curr_dist[w]){
              next_dist[v] = min(curr_dist[w] + 1, next_dist[v]);
            }
          }
          que.push(make_pair(next_dist[v], v));
        }
      }
    
      // Step 3.b
      while (!que.empty()){
        int d = que.top().first;
        int v = que.top().second; que.pop();
        if (d > next_dist[v]) continue;
        for (int w : fadj->at(v)){
          if (next_dist[w] > curr_dist[w] && next_dist[v] + 1 < next_dist[w]){
            next_dist[w] = next_dist[v] + 1;
            que.push(make_pair(next_dist[w], w));
          }
        }
      }

      for (int v : upd_nodes){
        CHECK(next_dist[v] != -1);
        if (next_dist[v] > curr_dist[v]) chg_nodes.push_back(v);
        curr_dist[v] = min(INF, next_dist[v]);
        next_dist[v] = -1;
      }
    }

    ReachabilityQuerier::ReachabilityQuerier(int source,
                                             int target,
                                             vector<vector<int> > *fadj,
                                             vector<vector<int> > *badj,
                                             SpecialPurposeReachabilityIndex *spr_index)
      : source(source), target(target), fadj(fadj), badj(badj), spr_index(spr_index)
    {
      distance.set_empty_key(-1);
      distance.set_deleted_key(-2);
      Build();
      source_in_mask  = spr_index->GetInMask(source);
      source_out_mask = spr_index->GetOutMask(source);
      target_in_mask  = spr_index->GetInMask(target);
      target_out_mask = spr_index->GetOutMask(target);
    }

    void ReachabilityQuerier::Build(){
      distance.clear();
      if (!ReachByTrees()){
        queue<int> que;
        que.push(source);
        distance[source] = 0;
        while (!que.empty()){
          int v = que.front(); que.pop();
          int d = distance[v];
          for (int w : fadj->at(v)){
            if (!this->Prune(w) && distance.find(w) == distance.end()){
              distance[w] = d + 1;
              que.push(w);
            }
          }
        }
      }
      UpdateMask();
    }

    void ReachabilityQuerier::InsertEdge(int u, int v){
      bool prev_tree_reach = (bool)(source_in_mask & target_out_mask);
      bool mask_change     = TargetMaskChanged();
      UpdateMask();
    
      if (ReachByTrees()){
        distance.clear();
      } else if (mask_change || prev_tree_reach){
        assert(!prev_tree_reach);
        this->Build();
      } else if (distance.find(u) != distance.end()){
        queue<pair<int, int> > que;
        if (GetDistance(v) > GetDistance(u) + 1){
          distance[v] = GetDistance(u) + 1;
          que.push(make_pair(v, distance[v]));
          // change_vs_ei++;
        }
        while (!que.empty()){
          int v = que.front().first;
          int d = que.front().second; que.pop();
          for (int w : fadj->at(v)){
            if (!this->Prune(w) && GetDistance(w) > d + 1){
              distance[w] = d + 1;
              que.push(make_pair(w, d + 1));
              // change_vs_ei++;
            }
          }
        }
      }
    }

    int ReachabilityQuerier::FindParent(int v) const {
      const auto &temp_dist = spr_index->temp_array;

      int dv = GetDistance(v);
    
      for (int w : badj->at(v)){
        bool not_red = temp_dist.at(w) == -1 || temp_dist.at(w) == GetDistance(w);
        if (not_red && dv == GetDistance(w) + 1){
          return w;
        }
      }
      return -1;
    }

    void ReachabilityQuerier::CollectChanges(const vector<int> &start_nodes, vector<int> &upd_nodes){
      auto &temp_dist = spr_index->temp_array;
    
      queue<int>  que;
      for (int v : start_nodes){
        if (FindParent(v) == -1){
          que.push(v);
          temp_dist.at(v) = INF;
          upd_nodes.push_back(v);
        }
      }
    
      while (!que.empty()){
        int v  = que.front(); que.pop();
        int dv = GetDistance(v);
        for (int w : fadj->at(v)){
          if ( temp_dist.at(w) != -1 || GetDistance(w) != dv + 1) continue;
        
          if (FindParent(w) == -1){
            temp_dist.at(w) = INF;
            que.push(w);
          } else {
            temp_dist.at(w) = GetDistance(w);
          }
          upd_nodes.push_back(w);
        }
      }
    }

    void ReachabilityQuerier::FixChanges(const vector<int> &upd_nodes){
    typedef pair<int, int> PI;
    priority_queue<PI, vector<PI>, greater<PI> > que;

    auto &temp_dist = spr_index->temp_array;
    
    for (int v : upd_nodes){
      if (temp_dist.at(v) <= GetDistance(v)) continue;
            
      for (int w : badj->at(v)){
        int dw = GetDistance(w);
        if (temp_dist.at(w) == -1 || temp_dist.at(w) == dw){
          temp_dist.at(v) = min(dw + 1, temp_dist.at(v));
        }
      }
      que.push(make_pair(temp_dist.at(v), v));
    }
        
    while (!que.empty()){
      int d = que.top().first;
      int v = que.top().second; que.pop();
      if (d > temp_dist.at(v)) continue;
        
      for (int w : fadj->at(v)){
        int dw = temp_dist.at(w);
        if (dw != -1 && dw > GetDistance(w) && d + 1 < dw){
          temp_dist.at(w) = temp_dist.at(v) + 1;
          que.push(make_pair(temp_dist.at(w), w));
        }
      }
    }
    
    for (int v : upd_nodes){
      if (temp_dist.at(v) < INF){
        distance[v] = temp_dist.at(v);
      } else {
        distance.erase(v);
      }
      temp_dist.at(v) = -1;
    }
    distance.resize(0);
  }

  void ReachabilityQuerier::DeleteEdge(int u, int v){
    bool prev_tree_reach = (bool)(source_in_mask & target_out_mask);
    bool mask_change = TargetMaskChanged();
    UpdateMask();
    
    if (ReachByTrees()){
      distance.clear();
    } else if (mask_change || prev_tree_reach){
      this->Build();
    } else {
      if (GetDistance(u) + 1 == GetDistance(v)){
        vector<int> upd_nodes;
        CollectChanges({v}, upd_nodes);
        FixChanges(upd_nodes);
        // change_vs_ed += upd_nodes.size();
        upd_nodes.clear();
      }
      
      typedef pair<int, int> PI;
      priority_queue<PI, vector<PI>, greater<PI> > que;
      const vector<int> *chg_nodes = spr_index->GetRCNodes();
      // change_vs_ed += chg_nodes->size();
      
      for (auto it = chg_nodes->begin(); it != chg_nodes->end(); it++){
        int v = *it;
        if (!this->Prune(v)){
          int best_dist = INF;
          for (int w : badj->at(v)){
            best_dist = min(best_dist, this->GetDistance(w));
          }
          
          if (best_dist < INF){
            distance[v] = best_dist + 1;
            que.push(PI(distance[v], v));
          }
        }
      }
      while (!que.empty()){
        int v = que.top().second; 
        int d = que.top().first; que.pop();
        if (d > distance[v]) continue;

        for (int w : fadj->at(v)){
          if (!this->Prune(w) && d + 1 < this->GetDistance(w)){
            distance[w] = d + 1;
            que.push(PI(distance[w], w));
          }
        }
      }
    }
  }

  void ReachabilityQuerier::InsertNode(int u){
    // 何もしない
    assert(0 <= u && (size_t)u < fadj->size());
    assert(0 <= u && (size_t)u < badj->size());
  }

  void ReachabilityQuerier::DeleteNode(int u, const vector<int> &u_out, const vector<int> &){
    assert(fadj->at(u).empty());
    assert(badj->at(u).empty());
    bool prev_tree_reach = (bool)(source_in_mask & target_out_mask);
    bool mask_change     = TargetMaskChanged();
    
    if (ReachByTrees()){
      distance.clear();
    } else if (mask_change || prev_tree_reach){
      this->Build();
    } else {
      if (distance.count(u) > 0){
        vector<int> start_nodes;
        vector<int> upd_nodes;
        int dist_u = GetDistance(u);
        for (int v : u_out){
          if (dist_u + 1 == GetDistance(v)) start_nodes.push_back(v);
        }
        distance.erase(u);
        CollectChanges(start_nodes, upd_nodes);
        FixChanges(upd_nodes);
        // change_vs_vd += upd_nodes.size();
        assert(distance.count(u) == 0);
      }
      
      typedef pair<int, int> PI;
      priority_queue<PI, vector<PI>, greater<PI> > que;
      const vector<int> *chg_nodes = spr_index->GetRCNodes();
      // change_vs_vd += chg_nodes->size();
      
      for (auto it = chg_nodes->begin(); it != chg_nodes->end(); it++){
        int v = *it;
        if (!this->Prune(v)){
          int best_dist = INF;
          for (int w : badj->at(v)){
            best_dist = min(best_dist, this->GetDistance(w));
          }
          
          if (best_dist < INF){
            distance[v] = best_dist + 1;
            que.push(PI(distance[v], v));
          }
        }
      }
      while (!que.empty()){
        int v = que.top().second; 
        int d = que.top().first; que.pop();
        if (d > distance[v]) continue;

        for (int w : fadj->at(v)){
          if (!this->Prune(w) && d + 1 < this->GetDistance(w)){
            distance[w] = d + 1;
            que.push(PI(distance[w], w));
          }
        }
      }
    }
  }

  bool ReachabilityQuerier::TargetMaskChanged() const {
    if (target_in_mask != spr_index->GetInMask(target)) return true;
    if (target_out_mask != spr_index->GetOutMask(target)) return true;
    return false;
  }

  bool ReachabilityQuerier::SourceMaskChanged() const {
    if (source_in_mask != spr_index->GetInMask(source)) return true;
    if (source_out_mask != spr_index->GetOutMask(source)) return true;
    return false;
  }
  
  void ReachabilityQuerier::UpdateMask(){
    source_in_mask  = spr_index->GetInMask(source);
    source_out_mask = spr_index->GetOutMask(source);
    target_in_mask  = spr_index->GetInMask(target);
    target_out_mask = spr_index->GetOutMask(target);
  }
  
  bool ReachabilityQuerier::Reach() const {
    return ReachByTrees() || distance.find(target) != distance.end();
  }

  const vector<int> ReachabilityQuerier::GetIndexNodes() const {
    vector<int> res;
    for (auto p : distance) res.push_back(p.first);
    return res;
  }

    
  
    SpecialPurposeReachabilityIndex::SpecialPurposeReachabilityIndex(vector<vector<int> >  *fadj, vector<vector<int> >  *badj, int num_rs)
      : fadj(fadj), badj(badj), id_manager(fadj->size()), num_rs(num_rs)
    {
      CHECK(fadj != nullptr && badj != nullptr && num_rs <= num_rs_limit);
    
      V = fadj->size();
      // #ifdef NDEBUG
      // JLOG_ADD_BENCHMARK("spr_index.construct_time"){
      // #endif
      has_change.resize(V, false);
      temp_array.resize(V, -1);
      for (int i = 0; i < 2; i++){
        reach_mask[i].resize(V, 0);
      }

      for (int k = 0; k < num_rs; k++){
        int root = rand() % V;
        roots.push_back(root);
        spts[0].push_back(new DynamicSPT(root, fadj, badj));
        spts[1].push_back(new DynamicSPT(root, badj, fadj));

        for (int v = 0; v < V; v++){
          for (int i = 0; i < 2; i++){
            SetBit(reach_mask[i][v], k, spts[i].back()->GetDistance(v) < INF);
          }
        }
      }
      // #ifdef NDEBUG
      // }
      // #endif 
    }
  
    SpecialPurposeReachabilityIndex::~SpecialPurposeReachabilityIndex(){
      for (auto prq : pr_queriers){
        delete prq;
      }
      
      for (int i = 0; i < 2; i++){
        for (auto spt : spts[i]) delete spt;
      }
    }
  
    void SpecialPurposeReachabilityIndex::InsertEdge(int u, int v){
      chg_nodes.clear();

      for (int k = 0; k < num_rs; k++){
        spts[0][k]->InsertEdge(u, v);
        spts[1][k]->InsertEdge(v, u);
      }
      CollectRCNodes();
      for (auto prq : pr_queriers){
        prq->InsertEdge(u, v);
      }
    }
  
    void SpecialPurposeReachabilityIndex::DeleteEdge(int u, int v){
      chg_nodes.clear();
      for (int k = 0; k < num_rs; k++){
        spts[0][k]->DeleteEdge(u, v);
        spts[1][k]->DeleteEdge(v, u);
      }
      CollectRCNodes();

      for (auto prq : pr_queriers){
        prq->DeleteEdge(u, v);
      }
    }
    
    void SpecialPurposeReachabilityIndex::InsertNode(int u){
      int new_V = max(u + 1, V);
      
      // Resize variables
      if (new_V > V){
        for (; V < new_V; V++){
          has_change.push_back(false);
          temp_array.push_back(-1);

          for (int i = 0; i < 2; i++){
            reach_mask[i].push_back(0);
          }
        }
        chg_nodes.clear();
        for (int k = 0; k < num_rs; k++){
          for (int i = 0; i < 2; i++){
            spts[i][k]->InsertNode(u);
          }
        }
        CollectRCNodes();
      }

      id_manager.MakeAlive(u);
      for (auto prq : pr_queriers){
        CHECK(prq != nullptr);
        prq->InsertNode(u);
      }
    }
  
    void SpecialPurposeReachabilityIndex::DeleteNode(int u, const vector<int> &u_out, const vector<int> &u_in){
      CHECK(fadj->at(u).empty() && badj->at(u).empty());
      chg_nodes.clear();
      id_manager.MakeDead(u);
      for (int k = 0; k < num_rs; k++){
        if (roots[k] == u){
          while (roots[k] == u && id_manager.NumAlive()) {
            roots[k] = id_manager.SampleAlive();
          }
          spts[0][k]->ChangeRoot(roots[k]);
          spts[1][k]->ChangeRoot(roots[k]);

          for (int v = 0; v < V; v++){
            for (int i = 0; i < 2; i++){
              SetBit(reach_mask[i][v], k, spts[i][k]->GetDistance(v) < INF);
            }
          }
        } else {
          spts[0][k]->DeleteNode(u, u_out, u_in);
          spts[1][k]->DeleteNode(u, u_in, u_out);
        }
      }
      CollectRCNodes();

      for (auto prq : pr_queriers){
        CHECK(prq != nullptr);
        prq->DeleteNode(u, u_out, u_in);
      }
    }

    ReachabilityQuerier *SpecialPurposeReachabilityIndex::CreateQuerier(int source, int target){
      ReachabilityQuerier *prq = new ReachabilityQuerier(source, target, fadj, badj, this);
      pr_queriers.push_back(prq);
      return prq;
    }
  
    const vector<pair<int, vector<int> > > SpecialPurposeReachabilityIndex::GetTrees() const {
      vector<pair<int, vector<int> > >  res;

      for (int k = 0; k < num_rs; k++){
        res.push_back(make_pair( roots[k], *spts[0][k]->GetTreeNodes()));
        res.push_back(make_pair(-roots[k], *spts[1][k]->GetTreeNodes()));
      }
      return res;
    }

    void SpecialPurposeReachabilityIndex::CollectRCNodes() {
      CHECK(chg_nodes.empty());

      for (int k = 0; k < num_rs; k++){
        for (int i = 0; i < 2; i++){
          const vector<int> *upd_nodes = spts[i][k]->GetDCNodes();
          for (auto it = upd_nodes->begin(); it != upd_nodes->end(); it++){
            int  v         = *it;
            bool cur_reach = GetBit(reach_mask[i][v], k);
            bool nxt_reach = spts[i][k]->GetDistance(v) < INF;
            if (cur_reach != nxt_reach && !has_change[v]) {
              chg_nodes.push_back(v);
              has_change[v] = true;
            }
            SetBit(reach_mask[i][v], k, nxt_reach);
          }
        }
      }
      
      for (int v : chg_nodes){
        has_change[v] = false;
      }
    }
  } /* special_purpose_reachability_index */
} /* dynamic_centrality */





