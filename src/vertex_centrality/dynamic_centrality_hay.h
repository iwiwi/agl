#ifndef DYNAMIC_CENTRALITY_HAY_H
#define DYNAMIC_CENTRALITY_HAY_H

#include "common.h"
#include "dynamic_centrality.h"
#include "hyper_edge.h"
#include "special_purpose_reachability_index.h"
#include <vector>
#include <cstdlib>
#include <queue>
using std::vector;
using std::pair;

namespace betweenness_centrality {
  
  class DynamicCentralityHAY : public DynamicCentralityBase {
  private:
    bool debug_mode;
    int num_samples;
    int tradeoff_param;
    vector<double>     score;
    vector<HyperEdge*> hyper_edges;

    // maintain ids that are assigned to each vertex.
    IDManager *id_manager;
    
    // temporal variables for index construction
    vector<int>    tmp_dist[2];
    vector<double> tmp_count[2];
    vector<int>    tmp_passable;
    
    // keep disjoint set union of nodes
    special_purpose_reachability_index::SpecialPurposeReachabilityIndex *spr_index;
    
    void Init();  // Initialize the arrays
    void Clear(); // Delete the array
    int SampleVertex() const ;
    vector<pair<int, int> > SampleVertexPairs() const ;
    
    bool InsertEdgeIntoGraph(int s, int t);
    bool DeleteEdgeFromGraph(int s, int t);
    bool InsertNodeIntoGraph(int v);
    bool DeleteNodeFromGraph(int v);
    inline bool ValidNode(int v) const { return vertex2id.count(v); }
    
  public:
    DynamicCentralityHAY() : debug_mode(false), tradeoff_param(0), id_manager(nullptr), spr_index(nullptr) { }
    ~DynamicCentralityHAY(){ Clear(); }
    
    virtual void PreCompute(const vector<pair<int, int> > &es, int num_samples);
    
    virtual double QueryCentrality(int v) const {
      size_t num_vs = vertex2id.size();
      return ValidNode(v) ? score[vertex2id.at(v)] / hyper_edges.size() * num_vs * num_vs : 0.0;
    }
    
    virtual void InsertEdge(int s, int t);
    virtual void DeleteEdge(int s, int t);
    virtual void InsertNode(int v);
    virtual void DeleteNode(int v);
    
    void SetTradeOffParam(int x) { tradeoff_param = x;}
    friend class HyperEdge;
  };
};

#endif /* DYNAMIC_CENTRALITY_HAY_H */


