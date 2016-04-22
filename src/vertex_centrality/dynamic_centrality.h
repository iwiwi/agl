#ifndef DYNAMIC_CENTRALITY_H
#define DYNAMIC_CENTRALITY_H

#include "centrality.h"

namespace betweenness_centrality {
  class DynamicCentralityBase
    : public betweenness_centrality::CentralityBase { 
  public:
    virtual double QueryCentrality(int v) const = 0;
    virtual void InsertNode(int v) = 0;
    virtual void DeleteNode(int v) = 0;
    virtual void InsertEdge(int u, int v) = 0;
    virtual void DeleteEdge(int u, int v) = 0;
  };
}

#endif /* DYNAMIC_CENTRALITY_H */
