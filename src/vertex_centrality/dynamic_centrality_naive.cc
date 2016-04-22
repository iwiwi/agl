#include "dynamic_centrality_naive.h"
#include "common.h"
using namespace std;

namespace betweenness_centrality {
  void DynamicCentralityNaive::PreCompute(const vector<pair<int, int> > &es_, int){
    es = es_;
    BuildGraph(es);
    for (auto &e : es){
      e.fst = vertex2id[e.fst];
      e.snd = vertex2id[e.snd];
    }
    cn.PreCompute(es);
  }

  void DynamicCentralityNaive::InsertNode(int v){
    CHECK(vertex2id.count(v) == 0);
    int V = vertex2id.size();
    vertex2id[v] = V;
  }

  void DynamicCentralityNaive::DeleteNode(int v){
    CHECK(vertex2id.count(v));
    int w = vertex2id[w];
    vertex2id.erase(v);
    
    size_t m = 0;
    for (size_t i = 0; i < es.size(); i++){
      if (es[i].fst != w && es[i].snd != w)  es[m++] = es[i];
    }
    
    if (m < es.size()){
      es.resize(m);
      cn.PreCompute(es);
    }
  }

  void DynamicCentralityNaive::InsertEdge(int u, int v){
    CHECK(vertex2id.count(u) && vertex2id.count(v));
    u = vertex2id[u];
    v = vertex2id[v];
    for (const auto &e : es){
      CHECK(e != make_pair(u, v));
    }
    es.emplace_back(u, v);
    cn.PreCompute(es);
  }

  void DynamicCentralityNaive::DeleteEdge(int u, int v){
    CHECK(vertex2id.count(u) && vertex2id.count(v));
    u = vertex2id[u];
    v = vertex2id[v];
    
    size_t m = 0;
    for (size_t i = 0; i < es.size(); i++){
      if (es[i] != make_pair(u, v)) es[m++] = es[i];
    }

    if (m < es.size()){
      es.resize(m);
      cn.PreCompute(es);
    }
  }
}
