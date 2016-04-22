#ifndef HYPEREDGE_H
#define HYPEREDGE_H

#include "common.h"
#include "special_purpose_reachability_index.h"
#include "sparsehash/dense_hash_map"
using std::vector;

namespace betweenness_centrality {

  template <typename T, typename E> using hash_map = google::dense_hash_map<T,E>;
  
  class DynamicCentralityHAY;
  
  class Ball {
  private:
    int source;
    int radius;
    vector<int> *tmp_dist;
    hash_map<int, int> distance;
    vector<vector<int> >  *fadj;
    vector<vector<int> >  *badj;
    
  public:
    Ball() : tmp_dist(nullptr) {
      distance.set_empty_key(-1);
      distance.set_deleted_key(-2);
    }
    void Build(const vector<std::pair<int, int> > &, vector<vector<int> > *, vector<vector<int> > *);
    void Trace(const vector<int> &start_nodes, vector<int> &dag_nodes);
    void DecreaseRadius();
    void InsertEdge(int u, int v);
    void DeleteEdge(int u, int v);
    void InsertNode(int v);  
    void DeleteNode(int u, const vector<int> &u_out, const vector<int> &u_in);
    void SetTempDist(vector<int> *tmp_dist){ this->tmp_dist = tmp_dist; }
    void UnsetTempDist(){ this->tmp_dist = nullptr; }
    
    inline bool HasNode(int v) const { return distance.find(v) != distance.end();}
    inline int GetDistance(int v) const {
      const auto iter = distance.find(v);
      assert(iter != distance.end());
      return iter->second;
    }
    inline int GetRadius() const { return radius; }
    inline size_t GetBallSize() const { return distance.size(); }
  private:  
    int FindParent(int v) const ;
    void CollectChanges(const vector<int> &start_nodes, vector<int> &upd_nodes);
    void FixChanges(const vector<int> &nodes);
  public:
    friend void Intersection(const Ball &, const Ball &, vector<int> &);
  };
  
  class HyperEdge {

  private: 
    bool is_connected;
    int  source;
    int  target;
    int  distance;
    Ball ball_s;
    Ball ball_t;
    hash_map<int, double> scores;
    hash_map<int, int>    dists;
    DynamicCentralityHAY *dch;
    special_purpose_reachability_index::ReachabilityQuerier *prq;
    
  public:
    HyperEdge(int s, int t, DynamicCentralityHAY *dch);
    ~HyperEdge(){ if (source != target && is_connected) SubWeight(); };
    void InsertEdge(int s, int t);
    void DeleteEdge(int u, int v);
    void InsertNode(int u);
    void DeleteNode(int u, const vector<int> &u_out, const vector<int> &u_in);
    
    inline int GetSource() const { return source; }
    inline int GetTarget() const { return target; }
    inline int ShortestPathLength() const { return distance; }
    inline bool IsConnected() const { return is_connected; }
    inline int GetNumNodes() const { return scores.size(); }

  private:
    bool BidirectionalSearch(int s, int t);
    void ComputeNumPaths(int s, const vector<vector<int> >  &adj, vector<int> &dist, vector<double> &count);
    bool RecomputeIndex();

    void UpdateDAGbyInsertion1(int u, int v);
    void UpdateDAGbyInsertion2(int u, int v);
    void UpdateDAGbyInsertion3(int u, int v);
    
    void CalcWeight();
    void CalcWeight(const vector<int> &dag_nodes);
    void AddWeight();
    void SubWeight();
  };
  

} /* betweenness_centrality */

#endif /* HYPEREDGE_H */
