#ifndef SPECIAL_PUPOSE_REACHABILITY_INDEX_H
#define SPECIAL_PUPOSE_REACHABILITY_INDEX_H

#include <vector>
#include <cstdlib>
#include <cassert>
#include <limits>
#include <numeric>
#include <iostream>
#include "sparsehash/dense_hash_map"
#include "id_manager.h"
using std::vector;

namespace betweenness_centrality {
  
  namespace special_purpose_reachability_index {
    class DynamicSPT;
    class ReachabilityQuerier;
    template <typename T, typename E> using hash_map = google::dense_hash_map<T,E>;
    constexpr static int INF = std::numeric_limits<int>::max() / 2;
    
    class SpecialPurposeReachabilityIndex {
    private:
      friend class ReachabilityQuerier;
      
      vector<int> roots;
      vector<int> reach_mask[2];
      vector<DynamicSPT*> spts[2];
      vector<vector<int> >  *fadj;
      vector<vector<int> >  *badj;
      IDManager id_manager;
    
      vector<int> temp_array;
      vector<int> has_change;
      vector<int> chg_nodes;
      vector<ReachabilityQuerier*> pr_queriers;
      int num_rs;
      int V;
      static const int num_rs_limit = 30;
    
    public:
    
      SpecialPurposeReachabilityIndex(vector<vector<int> >  *fadj, vector<vector<int> >  *badj, int num_rs);
      virtual ~SpecialPurposeReachabilityIndex();
      void InsertEdge(int u, int v);
      void DeleteEdge(int u, int v);
      void InsertNode(int u);
      void DeleteNode(int u, const vector<int> &u_out, const vector<int> &u_in);
      ReachabilityQuerier *CreateQuerier(int source, int target); 
      const vector<int> GetRoots() const { return roots; }
      const vector<std::pair<int, vector<int> > > GetTrees() const;
    
    private: 
      inline int GetOutMask(int v) const { assert(ValidNode(v)); return reach_mask[0][v]; }
      inline int GetInMask(int v) const { assert(ValidNode(v)); return reach_mask[1][v]; }
      inline bool ValidNode(int v) const { return 0 <= v && v < V; }
      const vector<int> *GetRCNodes() const { return &chg_nodes; }
      void CollectRCNodes();
    };


    class DynamicSPT {
    private: 
      int root;
      vector<int> curr_dist;
      vector<int> next_dist;
      vector<vector<int> >  *fadj;
      vector<vector<int> >  *badj;
      vector<int> chg_nodes;
    public:
      DynamicSPT(int r, vector<vector<int> >  *fadj, vector<vector<int> >  *badj);
      inline int  GetRoot() const { return root; }
      inline int  GetDistance(int v) { assert(ValidNode(v)); return curr_dist[v]; }
    
      void   ChangeRoot(int r);
      void   InsertEdge(int u, int v);
      void   DeleteEdge(int u, int v);
      void   InsertNode(int u);
      void   DeleteNode(int u, const vector<int> &u_out, const vector<int> &u_in);
      const vector<int> *GetDCNodes() { return &chg_nodes; }
      const vector<int> *GetTreeNodes() { return &curr_dist; }
      
    private:
      void Build();
      inline int  FindParent(int v);
      inline bool ValidNode(int v) { return 0 <= v && (size_t)v < fadj->size();}
      void CollectChanges(const vector<int> &start_nodes, vector<int> &upd_nodes);
      void FixChanges(const vector<int> &upd_nodes); 
    };

    class ReachabilityQuerier {
      friend class SpecialPurposeReachabilityIndex;
      int source;
      int target;
      int source_in_mask, source_out_mask;
      int target_in_mask, target_out_mask;
      hash_map<int, int>  distance;
      vector<vector<int> > *fadj;
      vector<vector<int> > *badj;
      SpecialPurposeReachabilityIndex   *spr_index;
      
    public: 
      ReachabilityQuerier(int source,
                          int tagret,
                          vector<vector<int> > *fadj,
                          vector<vector<int> > *badj,
                          SpecialPurposeReachabilityIndex *spr_index);
      bool Reach() const;
      inline int GetSource() const { return source; }
      inline int GetTarget() const { return target; }
      const vector<int> GetIndexNodes() const;
    private:
      void Build();
      void InsertEdge(int u, int v);
      void DeleteEdge(int u, int v);
      void InsertNode(int u);
      void DeleteNode(int u, const vector<int> &u_out, const vector<int> &u_in);
      bool TargetMaskChanged() const ;
      bool SourceMaskChanged() const ;
      void UpdateMask();
    
      inline bool ReachByTrees() const {
        return spr_index->GetOutMask(target) & spr_index->GetInMask(source);
      }
      
      inline bool Prune(int v) const {
        return ~spr_index->GetOutMask(target) & spr_index->GetOutMask(v);
      }
      
      inline int GetDistance(int v) const {
        auto iter = distance.find(v);
        return iter == distance.end() ? INF : iter->second;
      }

      int FindParent(int v) const ;
      void CollectChanges(const vector<int> &start_nodes, vector<int> &upd_nodes);
      void FixChanges(const vector<int> &upd_nodes);
    };

  } /* special_purpose_reachability_index */
} /* betweenness_centrality */




#endif /* SPECIAL_PUPOSE_REACHABILITY_INDEX_H */
