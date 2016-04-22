// #include "dynamic_index.hpp"
#include "dynamic_centrality_hay.h"
using namespace std;
using namespace betweenness_centrality::special_purpose_reachability_index;

namespace betweenness_centrality {
  
  vector<pair<int, int> > DynamicCentralityHAY::SampleVertexPairs() const {
    vector<pair<int, int> > res;
    if (debug_mode){
      for (int s = 0; (size_t)s < V; s++){
        for (int t = 0; (size_t)t < V; t++){
          res.emplace_back(s, t);
        }
      }
    } else {
      for (int i = 0; i < num_samples; i++){
        int s = id_manager->SampleAlive();
        int t = id_manager->SampleAlive();
        res.emplace_back(s, t);
      }
    }
    return res;
  }

  // グラフ以外の部分を初期化
  void DynamicCentralityHAY::Init(){
    spr_index  = new SpecialPurposeReachabilityIndex(&G[0], &G[1], 10);
    id_manager = new IDManager(V);
    score      = vector<double>(V, 0);
    for (int i = 0; i < 2; i++){
      tmp_dist[i]  = vector<int>(V, -1);
      tmp_count[i] = vector<double>(V, 0);
    }
    tmp_passable.resize(V, false);
  }
  
  void DynamicCentralityHAY::Clear(){
    SafeDelete(spr_index);
    SafeDelete(id_manager);
    
    for (auto &index : hyper_edges) SafeDelete(index);
    G[0].clear();
    G[1].clear();
    score.clear();
    hyper_edges.clear();
    for (int i = 0; i < 2; i++){
      tmp_dist[i].clear();
      tmp_count[i].clear();
    }
    tmp_passable.clear();
  }
  
  void DynamicCentralityHAY::
  PreCompute(const vector<pair<int, int> > &es, int num_samples_){    
    Clear();
    this->num_samples = num_samples_;
    if (num_samples == -1){
      debug_mode = true;
    } else {
      debug_mode = false;
    }
    
    BuildGraph(es);
    Init();
    
    auto vertex_pairs = SampleVertexPairs();
    for (const auto &vp : vertex_pairs){
      hyper_edges.push_back(new HyperEdge(vp.fst, vp.snd, this));
    }
  }

  // 辺 {s, t}がすでにあった場合は何もせずfalseをかえす
  bool DynamicCentralityHAY::InsertEdgeIntoGraph(int s, int t){
    auto &f_adj = G[0];
    auto &b_adj = G[1];
    
    auto fiter = lower_bound(f_adj[s].begin(), f_adj[s].end(), t);
    if (fiter != f_adj[s].end() && *fiter == t){
      return false; 
    } else {
      auto biter = lower_bound(b_adj[t].begin(), b_adj[t].end(), s);
      f_adj[s].insert(fiter, t);
      b_adj[t].insert(biter, s);
      return true;
    }
  }
  
  // 辺 {s, t}が存在しない場合何もしないでfalseを返す
  bool DynamicCentralityHAY::DeleteEdgeFromGraph(int s, int t){
    auto &f_adj = G[0];
    auto &b_adj = G[1];
    auto fiter = lower_bound(f_adj[s].begin(), f_adj[s].end(), t);
    if (fiter != f_adj[s].end() && *fiter == t){
      auto biter = lower_bound(b_adj[t].begin(), b_adj[t].end(), s);
      f_adj[s].erase(fiter);
      b_adj[t].erase(biter);
      return true;
    } else {
      return false;  
    } 
  }

  bool DynamicCentralityHAY::InsertNodeIntoGraph(int v){
    if (vertex2id.count(v)){
      return false;
    } else {
      CHECK(id_manager->Full());
      if (id_manager->Full()){
        id_manager->Add(false);
        V++;
      }
      CHECK(id_manager->Size() == V);
      vertex2id[v] = id_manager->SampleDead();
      CHECK(vertex2id.size() == V);
      id_manager->MakeAlive(vertex2id[v]);

      while (G[0].size() < id_manager->Size()){
        G[0].push_back(vector<int>());
        G[1].push_back(vector<int>());
        score.push_back(0);
        for (int i = 0; i < 2; i++){
          tmp_dist[i].push_back(-1);
          tmp_count[i].push_back(0);
        }
        tmp_passable.push_back(false);
      }
      CHECK(G[0].size() == V);
      return true;
    }
  }

  bool DynamicCentralityHAY::DeleteNodeFromGraph(int u){
    vector<int> u_out(G[0][u]);
    vector<int> u_in(G[1][u]);

    for (int v : u_out){
      CHECK(u != v);
      auto iter = lower_bound(G[1][v].begin(), G[1][v].end(), u);
      CHECK(iter != G[1][v].end() && *iter == u);
      G[1][v].erase(iter);
    }
    
    for (int v : u_in){
      CHECK(u != v);
      auto iter = lower_bound(G[0][v].begin(), G[0][v].end(), u);
      CHECK(iter != G[0][v].end() && *iter == u);
      G[0][v].erase(iter);
    }
    G[0][u].clear();
    G[1][u].clear();
    return true;
  }

  void DynamicCentralityHAY::InsertEdge(int s, int t){
    CHECK(vertex2id.count(s) && vertex2id.count(t));
    s = vertex2id[s];
    t = vertex2id[t];
    if (InsertEdgeIntoGraph(s, t)){
      spr_index->InsertEdge(s, t);
      for (auto e : hyper_edges){
        e->InsertEdge(s, t);
      }
    }
  }

  void DynamicCentralityHAY::DeleteEdge(int s, int t){
    CHECK(vertex2id.count(s) && vertex2id.count(t));
    s = vertex2id[s];
    t = vertex2id[t];
    if (DeleteEdgeFromGraph(s, t)){
      spr_index->DeleteEdge(s, t);
      for (auto &e : hyper_edges){
        e->DeleteEdge(s, t);
      }
    }
  }
  
  void DynamicCentralityHAY::InsertNode(int u){
    if (InsertNodeIntoGraph(u) && ValidNode(u)){
      u = vertex2id[u];
      spr_index->InsertNode(u);
      CHECK(vertex2id.size() == V);
      if (debug_mode){
        for (int s = 0; size_t(s) < V; s++)
          for (int t = 0; size_t(t) < V; t++)
            if ((s == u && ValidNode(t)) || (ValidNode(s) && t == u))
              hyper_edges.push_back(new HyperEdge(s, t, this));
      } else {
        size_t n     = vertex2id.size();
        CHECK(n > 0);
        double prob1 = (double)(n - 1) / n * (n - 1) / n;;
        double prob2 = 1.0 / n / n;
        
        for (auto &e : hyper_edges){
          int new_source = -1;
          int new_target = -1;
          double q = (double)rand() / RAND_MAX;
          if (q < prob1){ // no change
            continue;
          } else if (q < prob2 + prob1){
            new_source = new_target = u; // re-sample pair of same vertices (really low prob.)
          } else {
            new_source = new_target = u;
            CHECK(n > 1u);
            while (new_target == u){
              new_target = id_manager->SampleAlive();
            }
            
            if (((double)rand() / RAND_MAX) < 0.5){
              swap(new_source, new_target);
            }
          }
          SafeDelete(e);
          e = new HyperEdge(new_source, new_target, this);
        }
      }
      for (auto e : hyper_edges){
        e->InsertNode(u);
      }
    }
  }
  
  void DynamicCentralityHAY::DeleteNode(int u){
    if (vertex2id.count(u) == 0){
      return;
    }
    
    int v = vertex2id[u];
    CHECK(id_manager->MakeDead(v));
    vertex2id.erase(u);
    CHECK(vertex2id.size() == id_manager->NumAlive());
    
    vector<int> v_out(G[0].at(v));
    vector<int> v_in (G[1].at(v));
    
    if (DeleteNodeFromGraph(v)){
      spr_index->DeleteNode(v, v_out, v_in);
      if (debug_mode){
        vector<HyperEdge*> new_hs;
        for (auto &e : hyper_edges) {
          if (e->GetSource() == v || e->GetTarget() == v) {
            SafeDelete(e);
          } else {
            e->DeleteNode(v, v_out, v_in);
            new_hs.push_back(e);
          }
        }
        hyper_edges = new_hs;
        
      } else {
        for (auto &e : hyper_edges){
          if (e->GetSource() == v || e->GetTarget() == v){
            int new_source = id_manager->SampleAlive();
            int new_target = id_manager->SampleAlive();
            CHECK(new_target != v && new_source != v);
            SafeDelete(e);
            e = new HyperEdge(new_source, new_target, this);
          } else {
            e->DeleteNode(v, v_out, v_in);
          }
        }
      }
      CHECK(score[v] < 1e-9);
    }
  }

} /* betweenness_centrality */
