#include "hyper_edge.h"
#include "dynamic_centrality_hay.h"
#include "common.h"
#include <set>
#include <queue>
using namespace std;

namespace betweenness_centrality {

  void Ball::Build(const vector<pair<int, int> > &nodes, vector<vector<int> > *fadj, vector<vector<int> > *badj){
    radius = 0;
    distance.clear();
    this->fadj = fadj;
    this->badj = badj;
    distance.insert(nodes.begin(), nodes.end());
    
    for (auto p : nodes){
      if (p.second == 0){
        source = p.first;
      }
      radius = max(radius, p.second);
    }
  }

  void Ball::Trace(const vector<int> &start_nodes, vector<int> &dag_nodes){
    // start_nodesがボールの中心から等しい距離にあることを仮定
    queue<int> que;
    set<int>   visit;
    for (int v : start_nodes){
      que.push(v);
    }
    while (!que.empty()){
      int v = que.front(); que.pop();
      int d = GetDistance(v);
      for (int w : badj->at(v)){
        if (HasNode(w) && GetDistance(w) == d - 1 && visit.count(w) == 0){
          visit.insert(w);
          dag_nodes.push_back(w);
          que.push(w);
        }
      }
    }
  }
  
  void Ball::DecreaseRadius(){
    if (radius > 0){
      for (auto iter = distance.begin(); iter != distance.end();){
        if (iter->second == radius){
          distance.erase(iter++);
        } else {
          iter++;
        }
      }
      distance.resize(0);
      radius--;
    }
  }
  
  void Ball::InsertEdge(int u, int v){
    if (HasNode(u) && GetDistance(u) < radius){
      queue<pair<int, int> > que;
      if (!HasNode(v) || GetDistance(v) > GetDistance(u) + 1){
        
        distance[v] = GetDistance(u) + 1;
        que.push(make_pair(v, distance[v]));
      }

      while (!que.empty()){
        int v = que.front().first;
        int d = que.front().second; que.pop();
        if (d == radius) continue;
        for (int w : fadj->at(v)){
          if (!HasNode(w) || GetDistance(w) > d + 1){
            distance[w] = d + 1;
            que.push(make_pair(w, d + 1));
          }
        }
      }
    }
  }

  int Ball::FindParent(int v) const {
    assert(HasNode(v));
    int dv = GetDistance(v);
    for (int w : badj->at(v)){
      if (!HasNode(w)) continue;
      int prev_d = GetDistance(w);
      int curr_d = tmp_dist->at(w) == -1 ? prev_d : tmp_dist->at(w);
      if (prev_d == curr_d && dv == curr_d + 1){
        return w;
      }
    }
    return -1;
  }
  
  void Ball::CollectChanges(const vector<int> &start_nodes, vector<int> &updated_nodes){
    queue<int>  que;
    const int INF = fadj->size() + 1;
        
    for (int v : start_nodes){
      assert(HasNode(v));
      if (FindParent(v) == -1){
        que.push(v);
        tmp_dist->at(v) = INF;
        updated_nodes.push_back(v);
      }
    }
        
    while (!que.empty()){
      int v  = que.front(); que.pop();
      int dv = GetDistance(v);
      if (dv >= radius) continue;
      
      for (int w : fadj->at(v)){
        assert(this->HasNode(w));
        if (tmp_dist->at(w) != -1 || GetDistance(w) != dv + 1) continue;
        
        if (FindParent(w) == -1){
          tmp_dist->at(w) = INF;
          que.push(w);
        } else {
          tmp_dist->at(w) = GetDistance(w);
        }
        assert(tmp_dist->at(w) != -1);
        updated_nodes.push_back(w);
      }
    }
  }

  void Ball::FixChanges(const vector<int> &nodes){
    // 辺・頂点削除後のnodesの中心からの距離を確定する
    typedef pair<int, int> PI;
    priority_queue<PI, vector<PI>, greater<PI> > que;
        
    for (int v : nodes){
      assert(tmp_dist->at(v) != -1);
      if (tmp_dist->at(v) <= GetDistance(v)) continue;
          
      for (int w : badj->at(v)){
        if (!this->HasNode(w)) continue;
        
        int dw = GetDistance(w);
        if (tmp_dist->at(w) == -1 || tmp_dist->at(w) == dw){
          tmp_dist->at(v) = min(dw + 1, tmp_dist->at(v));
        }
      }
      que.push(PI(tmp_dist->at(v), v));
    }
    
    while (!que.empty()){
      int d = que.top().first;
      int v = que.top().second; que.pop();
      if (d > tmp_dist->at(v) || d >= radius) continue;
      
      for (int w : fadj->at(v)){
        assert(this->HasNode(w));
        bool is_red = tmp_dist->at(w) != -1 && tmp_dist->at(w) > GetDistance(w);
        
        if (is_red && tmp_dist->at(v) + 1 < tmp_dist->at(w)){
          tmp_dist->at(w) = tmp_dist->at(v) + 1;
          que.push(PI(tmp_dist->at(w), w));
        }
      }
    }
    for (int v : nodes){
      if (tmp_dist->at(v) <= radius){
        distance[v] = tmp_dist->at(v);
      } else {
        distance.erase(v);
      }
      tmp_dist->at(v) = -1;
    }
    distance.resize(0);
  }
  
  void Ball::DeleteEdge(int u, int v){
    if (HasNode(u) && HasNode(v) && GetDistance(u) + 1 == GetDistance(v)){
      vector<int> start_nodes = {v};
      vector<int> upd_nodes;
      CollectChanges(start_nodes, upd_nodes);
      FixChanges(upd_nodes);
    }
  }
  
  void Ball::InsertNode(int v){
    assert(v < (int)fadj->size() && v < (int)badj->size());
    assert(tmp_dist == nullptr || v < (int)tmp_dist->size());
  }

  void Ball::DeleteNode(int u, const vector<int> &u_out, const vector<int> &){
    // 辺の情報の更新はボールの更新より後
    if (HasNode(u)){
      distance.erase(u);
      vector<int> start_nodes;
      vector<int> upd_nodes;

      for (int v : u_out){
        if (this->HasNode(v)) start_nodes.push_back(v);
      }
      CollectChanges(start_nodes, upd_nodes);
      FixChanges(upd_nodes);
    }
  }

  
  void Intersection(const Ball &bp, const Ball &bq, vector<int> &vec){
    vec.clear();
    const auto &bp_dist = bp.distance;
    const auto &bq_dist = bq.distance;
    for (const auto &p : bp_dist){
      if (bq_dist.find(p.first) != bq_dist.end()) vec.push_back(p.first);
    }
  }

  
  bool HyperEdge::
  BidirectionalSearch(int s, int t){
    assert(s != t);
    int         s_curr = 0, s_next = 2;
    int         t_curr = 1, t_next = 3;
    queue<int>  que[4];
    vector<int> update[2];
    que[s_curr].push(s); dch->tmp_dist[s_curr][s] = 0; update[s_curr].push_back(s);
    que[t_curr].push(t); dch->tmp_dist[t_curr][t] = 0; update[t_curr].push_back(t);
    
    bool found = false;
    while (!que[s_curr].empty() && !que[t_curr].empty()){
      int &curr = (update[0].size() <= update[1].size()) ? s_curr : t_curr;
      int &next = (update[0].size() <= update[1].size()) ? s_next : t_next;
      bool from_s = curr % 2 == 0;
      
      while (!que[curr].empty()){
        int v = que[curr].front(); que[curr].pop();
        int p = curr % 2;
        const auto &adj = from_s ? dch->G[0][v] : dch->G[1][v];
        
        for (int w : adj){
          int &src_d = dch->tmp_dist[    p][w];
          int &dst_d = dch->tmp_dist[1 - p][w];
          if (src_d != -1) continue;
          if (dst_d != -1) found = true;
          que[next].push(w);
          update[p].push_back(w);
          dch->tmp_dist[p][w] = dch->tmp_dist[p][v] + 1;
        }
      }
      if (found) goto LOOP_END;
      swap(curr, next);
    }
  LOOP_END:
    
    if (found){
      for (int i = 0; i < 2; i++){
        vector<pair<int, int> > nodes;
        for (auto v : update[i]){
          nodes.push_back(make_pair(v, dch->tmp_dist[i][v]));
        }
        if (i == 0){
          ball_s.Build(nodes, &dch->G[0], &dch->G[1]);
        } else {
          ball_t.Build(nodes, &dch->G[1], &dch->G[0]);
        }
      }
    } 

    for (int i = 0; i < 2; i++){
      for (int v : update[i]) dch->tmp_dist[i][v] =  -1;
    }
    return found;
  }

  void HyperEdge::
  ComputeNumPaths(int s, const vector<vector<int> >  &adj, vector<int> &dist, vector<double> &count){
    queue<int> que;
    que.push(s);
    dist[s] = 0;
    count[s] = 1;
    
    while (!que.empty()){
      int v = que.front(); que.pop();
      for (int w : adj[v]){
        int next_dist = dist[v] + 1;
        if (!dch->tmp_passable[w]) continue;
        if (dist[w] == -1){
          dist[w] = next_dist;
          que.push(w);
        }
        if (dist[w] == next_dist){
          count[w] += count[v];
        }
      }
    }
  }

  void HyperEdge::CalcWeight(){
    vector<int> common_nodes;
    Intersection(ball_s, ball_t, common_nodes);
    vector<int> dag_nodes = common_nodes;
    ball_s.Trace(common_nodes, dag_nodes);
    ball_t.Trace(common_nodes, dag_nodes);
    CalcWeight(dag_nodes);
  }
  
  void HyperEdge::CalcWeight(const vector<int> &dag_nodes){
    assert(is_connected);
    scores.clear();
    dists.clear();
    vector<int>    &dist_s  = dch->tmp_dist[0];
    vector<int>    &dist_t  = dch->tmp_dist[1];
    vector<double> &count_s = dch->tmp_count[0];
    vector<double> &count_t = dch->tmp_count[1];
    
    for (int v : dag_nodes){
      dch->tmp_passable[v] = true;
    }
    
    ComputeNumPaths(source, dch->G[0], dist_s, count_s);
    ComputeNumPaths(target, dch->G[1], dist_t, count_t);
    
    double num_paths = count_s[target];
    assert(num_paths > 0);

    this->distance = dist_s[target];
      
    for (int v : dag_nodes){
      double w  = count_s[v] * count_t[v] / num_paths;
      if (dist_s[v] + dist_t[v] == distance){
        this->scores[v] = w;
        this->dists[v]  = dist_s[v];
      }
      count_s[v] = count_t[v] = 0;
      dist_s[v] = dist_t[v] = -1;
      dch->tmp_passable[v] = false;
    }

    while (ball_s.GetRadius() + ball_t.GetRadius() + dch->tradeoff_param >= distance){
      auto &ball = ball_s.GetBallSize() > ball_t.GetBallSize() ? ball_s : ball_t;
      ball.DecreaseRadius();
      
      if (ball_s.GetRadius() + ball_t.GetRadius() == 0) break;
    }
  }
  
  void HyperEdge::AddWeight(){
    if (!is_connected) return;
    for (const auto p : scores){
      if (p.first != source && p.first != target){
        dch->score[p.first] += p.second;
      }
    }
  }
  
  void HyperEdge::SubWeight(){
    if (!is_connected) return;
    for (const auto p : scores){
      if (p.first != source && p.first != target){
        dch->score[p.first] -= p.second;
      }
    }
  }

  HyperEdge::HyperEdge(int s, int t, DynamicCentralityHAY *dch)
    : is_connected(false), source(s), target(t), dch(dch)
  {
    scores.set_empty_key(-1); scores.set_deleted_key(-2);
    dists.set_empty_key(-1); dists.set_deleted_key(-2);
    
    if (s != t){
      prq = dch->spr_index->CreateQuerier(s, t);
      is_connected = BidirectionalSearch(s, t);
      if (is_connected){
        CalcWeight();
        // cout << s << " " << t << " OK" << endl;
        // for (const auto &p : scores){
        //   cout << p << endl;
        // }
        AddWeight();
      }
    }
  }

  bool HyperEdge::RecomputeIndex(){
    scores.clear();
    dists.clear();
    is_connected = BidirectionalSearch(source, target);
    if (is_connected){
      CalcWeight();
    }
    return is_connected;
  }
  
  void Explore(int start_node,
               int max_radius,
               const Ball &goal,
               const vector<vector<int> >  &fadj,
               const vector<vector<int> >  &badj, 
               vector<int> &tmp_dist,
               vector<int> &dag_nodes,
               vector<int> &intersection)
  {
    if (goal.HasNode(start_node)){
      intersection.push_back(start_node);
      return;
    }
    
    if (max_radius == 0){
      return;
    }
    
    vector<int> updated_nodes;
    queue<int>  que;

    que.push(start_node);
    updated_nodes.push_back(start_node);
    tmp_dist[start_node] = 0;

    while (!que.empty()){
      int v = que.front(); que.pop();
      
      if (tmp_dist[v] == max_radius) break;
      
      for (int w : fadj[v]){
        if (tmp_dist[w] == -1){
          tmp_dist[w] = tmp_dist[v] + 1;
          que.push(w);
          updated_nodes.push_back(w);
          
          if (goal.HasNode(w)){
            max_radius = min(max_radius, tmp_dist[w]);
            intersection.push_back(w);
          }
        }
      }
    }

    if (!intersection.empty()){
      queue<int> que;
      for (int v: intersection){
        que.push(v);
      }
      set<int> visited_nodes;

      while (!que.empty()){
        int v = que.front(); que.pop();
        for (int w: badj[v]){
          if (visited_nodes.count(w) == 0 && tmp_dist[w] + 1 == tmp_dist[v] && tmp_dist[w] >= 0){
            que.push(w);
            dag_nodes.push_back(w);
            visited_nodes.insert(w);
          }
        }
      }
    }
    
    for (int v : updated_nodes){
      tmp_dist[v] = -1;
    }
  }

  

  void HyperEdge::UpdateDAGbyInsertion1(int u, int v){
    assert(ball_s.HasNode(u));

    bool u_in_s = ball_s.HasNode(u);
    bool v_in_s = ball_s.HasNode(v);

    int du = u_in_s ? ball_s.GetDistance(u) : dch->V + 1;
    int dv = v_in_s ? ball_s.GetDistance(v) : dch->V + 1;
    
    
    if (!v_in_s || du + 1 == dv){
          
      int max_radius = v_in_s ?
        (distance - dv - ball_t.GetRadius()) :
        (distance - du - ball_t.GetRadius());
      vector<int> &dist_s = dch->tmp_dist[0];
      vector<int>  dag_nodes;
      vector<int>  inter_nodes;
          
      Explore(v, max_radius, ball_t, dch->G[0], dch->G[1],
              dist_s, dag_nodes, inter_nodes);
          
      // compute DAG!
      if (!inter_nodes.empty()){
        ball_s.Trace({u}        , dag_nodes);
        ball_t.Trace(inter_nodes, dag_nodes);
        
        for (int w : inter_nodes){
          dag_nodes.push_back(w);
        }
        dag_nodes.push_back(u);
        
        for (auto p : scores){
          dag_nodes.push_back(p.first);
        }
        sort(dag_nodes.begin(), dag_nodes.end());
        dag_nodes.erase(unique(dag_nodes.begin(), dag_nodes.end()), dag_nodes.end());
        
        SubWeight();
        is_connected = RecomputeIndex();
        CalcWeight(dag_nodes);
        AddWeight();
      }
    }
  }
  

  void HyperEdge::UpdateDAGbyInsertion2(int u, int v){
    // // UpdateDAGbyInsertion1との重複がひどい
    assert(ball_t.HasNode(v));
    
    bool u_in_t = ball_t.HasNode(u);
    bool v_in_t = ball_t.HasNode(v);

    int du = u_in_t ? ball_t.GetDistance(u) : dch->V + 1;
    int dv = v_in_t ? ball_t.GetDistance(v) : dch->V + 1;
    
    if (!u_in_t || dv + 1 == du){
      
      int max_radius = u_in_t ?
        (distance - du - ball_s.GetRadius()) :
        (distance - dv - ball_s.GetRadius());
      vector<int> &dist_t = dch->tmp_dist[0];
      vector<int>  dag_nodes;
      vector<int>  inter_nodes;
      
      Explore(u, max_radius, ball_s, dch->G[1], dch->G[0],
              dist_t, dag_nodes, inter_nodes);
          
      // compute DAG!
      if (!inter_nodes.empty()){
        ball_t.Trace({v}        , dag_nodes);
        ball_s.Trace(inter_nodes, dag_nodes);
        for (int w : inter_nodes){
          dag_nodes.push_back(w);
        }
        dag_nodes.push_back(v);
        for (auto p : scores){
          dag_nodes.push_back(p.first);
        }
        sort(dag_nodes.begin(), dag_nodes.end());
        dag_nodes.erase(unique(dag_nodes.begin(), dag_nodes.end()), dag_nodes.end());
        
        SubWeight();
        CalcWeight(dag_nodes);
        AddWeight();
      }
    }
  }

  void HyperEdge::UpdateDAGbyInsertion3(int u, int v){
    // // UpdateDAGbyInsertion1との重複がひどい
    bool u_in_s = ball_s.HasNode(u), v_in_s = ball_s.HasNode(v);
    bool u_in_t = ball_t.HasNode(u), v_in_t = ball_t.HasNode(v);
    assert(!u_in_s && !u_in_t && !v_in_s && !v_in_t);
    
    vector<int> dag_nodes;
    vector<int> inter_nodes1;
    vector<int> inter_nodes2;
    const auto &fadj = dch->G[0];
    const auto &badj = dch->G[1];
    vector<int> &dist_s = dch->tmp_dist[0];
    vector<int> &dist_t = dch->tmp_dist[1];

    int max_radius = dch->tradeoff_param;
    Explore(u, max_radius, ball_s, badj, fadj, dist_t, dag_nodes, inter_nodes1);
    Explore(v, max_radius, ball_t, fadj, badj, dist_s, dag_nodes, inter_nodes2);
    
    if (!inter_nodes1.empty() && !inter_nodes2.empty()){
      // compute DAG!
      ball_s.Trace(inter_nodes1, dag_nodes);
      ball_t.Trace(inter_nodes2, dag_nodes);
      
      for (int w : inter_nodes1){
        dag_nodes.push_back(w);
      }
      
      for (int w : inter_nodes2){
        dag_nodes.push_back(w);
      }
      
      for (auto p : scores){
        dag_nodes.push_back(p.first);
      }

      sort(dag_nodes.begin(), dag_nodes.end());
      dag_nodes.erase(unique(dag_nodes.begin(), dag_nodes.end()), dag_nodes.end());

      SubWeight();
      CalcWeight(dag_nodes);
      AddWeight();
    }
  }
  
  void HyperEdge::InsertEdge(int u, int v){
    // // cout << "INSERT: " << u << " " << v << " " << source << " " << target << endl;
    if (source == target) return;
    
    if (is_connected){
      // 先にボールを更新する.
      if (ball_s.HasNode(u)){
        ball_s.SetTempDist(&dch->tmp_dist[0]);
        ball_s.InsertEdge(u, v);
        ball_s.UnsetTempDist();
      }
      
      if (ball_t.HasNode(v)){
        ball_t.SetTempDist(&dch->tmp_dist[0]);
        ball_t.InsertEdge(v, u);
        ball_t.UnsetTempDist();
      }
      
      bool u_in_s = ball_s.HasNode(u), v_in_s = ball_s.HasNode(v);
      bool u_in_t = ball_t.HasNode(u), v_in_t = ball_t.HasNode(v);
      
      if (u_in_s){
        UpdateDAGbyInsertion1(u, v);
      } else if (v_in_t){
        UpdateDAGbyInsertion2(u, v);
      } else if (!u_in_s && !v_in_s && !u_in_t && !v_in_t){
        UpdateDAGbyInsertion3(u, v);
      }
    } else if (prq->Reach()){
      // cout << "REACH: " << source << " " << target << endl;
      SubWeight();
      is_connected = RecomputeIndex();
      AddWeight();
    }
  }

  void HyperEdge::InsertNode(int u){
    // ここに来た時点でdch側はすでにグラフを更新している
    assert(0 <= u && size_t(u) < dch->G[0].size() && size_t(u) < dch->G[1].size());
    if (is_connected){
      ball_s.InsertNode(u);
      ball_t.InsertNode(u);
    }
  }

  void HyperEdge::
  DeleteEdge(int u, int v){
    if (source == target || !is_connected) return;
    
    auto u_iter = scores.find(u);
    auto v_iter = scores.find(v);
    bool dag_update = u_iter != scores.end() && v_iter != scores.end() && dists[u] + 1 == dists[v];
    
    if (dag_update){
      // DAGの更新が必要
      if (Equal(u_iter->second, 1.0) && Equal(v_iter->second, 1.0)){
        // sourceとtargetが非連結 <=> scores[u]=scores[v]=1
        SubWeight();
        is_connected = RecomputeIndex();
        AddWeight();
        // ボールの再計算も終わっている 
        return;
      } else {
        // DAG上の頂点集合が変わってしまうことに注意
        SubWeight();
        vector<int> dag_nodes;
        for (auto p : scores){
          dag_nodes.push_back(p.first);
        }
        CalcWeight(dag_nodes);
        AddWeight();
      }
    }
    
    ball_s.SetTempDist(&dch->tmp_dist[0]);
    ball_s.DeleteEdge(u, v);
    ball_s.UnsetTempDist();

    ball_t.SetTempDist(&dch->tmp_dist[1]);
    ball_t.DeleteEdge(v, u);
    ball_t.UnsetTempDist();
  }

  void HyperEdge::DeleteNode(int u, const vector<int> &u_out, const vector<int> &u_in){
    assert(dch->G[0][u].empty() && dch->G[1][u].empty());
    
    auto u_iter = scores.find(u);
    if (u_iter == scores.end()) return;
    
    if (Equal(u_iter->second, 1.0)){
      SubWeight();
      is_connected = RecomputeIndex();
      AddWeight();
      return;
    } else {
      // DAG上で再計算、ただしuは通らない
      SubWeight();
      vector<int> dag_nodes;
      for (auto p : scores){
        if (p.first != u) dag_nodes.push_back(p.first);
      }
      CalcWeight(dag_nodes);
      AddWeight();
    }
    ball_s.SetTempDist(&dch->tmp_dist[0]);
    ball_s.DeleteNode(u, u_out, u_in);
    ball_s.UnsetTempDist();
    
    ball_t.SetTempDist(&dch->tmp_dist[1]);
    ball_t.DeleteNode(u, u_in, u_out);
    ball_t.UnsetTempDist();
  }
  
  
} /* betweenness_centrality */

