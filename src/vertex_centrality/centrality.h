#ifndef CENTRALITY_H
#define CENTRALITY_H

#include "common.h"
#include <vector>
#include <unordered_map>
using std::vector;
using std::unordered_map;

namespace betweenness_centrality {

  class CentralityBase {
  protected:
    size_t V;
    size_t E;
    vector<vector<int> > G[2];
    unordered_map<int, int> vertex2id;
    void BuildGraph(const vector<std::pair<int, int> > &es);
    
  public:
    virtual ~CentralityBase(){};
    virtual void PreCompute(const vector<std::pair<int, int> > &es, int num_samples = -1) = 0;
    virtual double QueryCentrality(int v) const = 0;
  };
  
  class CentralityNaive : public CentralityBase {
    vector<double> centrality_map;
  public:
    virtual void PreCompute(const vector<std::pair<int, int> > &es, int num_samples = -1);
    virtual double QueryCentrality(int v) const {
      return vertex2id.count(v) ? centrality_map[vertex2id.at(v)] : 0;
    }
  };
  
  class CentralitySample : public CentralityBase {
    enum Direction {
      Forward,
      Backward
    };
    vector<double> centrality_map;
    vector<bool>   temp_on_DAG;
    vector<vector<int> > temp_distance;
    vector<vector<double> > temp_num_paths;
    
    vector<int> ComputeDAG(int source, int target);
    void BreadthFirstSearchOnDAG(int source, Direction dir);
  public:
    virtual void PreCompute(const vector<std::pair<int, int> > &es, int num_samples = -1);
    virtual double QueryCentrality(int v) const {
      return vertex2id.count(v) ? centrality_map[vertex2id.at(v)] : 0;
    }
  };
};



#endif /* CENTRALITY_H */
