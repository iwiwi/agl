#include "centrality.h"
#include <cassert>
#include <queue>
#include <algorithm>
using namespace std;

namespace betweenness_centrality {
  void CentralityBase::BuildGraph(const vector<pair<int, int> > &es){
    // CHECK(!es.empty());
    vertex2id.clear();
    G[0].clear();
    G[1].clear();
    
    if (!es.empty()){
      V = 0;
      E = es.size();
      for (const auto &e : es){
        if (vertex2id.count(e.fst) == 0) vertex2id[e.fst] = V++;
        if (vertex2id.count(e.snd) == 0) vertex2id[e.snd] = V++;
      }
      G[0].resize(V);
      G[1].resize(V);
      for (const auto &e : es){
        int u = vertex2id[e.fst];
        int v = vertex2id[e.snd];
        if (u != v){
          G[0][u].push_back(v);
          G[1][v].push_back(u);
        }
      }
    } else {
      V = 1;
      E = 0;
      G[0].resize(V);
      G[1].resize(V);
      vertex2id[0] = 0;
    }

    for (size_t v = 0; v < V; v++){
      sort(G[0][v].begin(), G[0][v].end());
      sort(G[1][v].begin(), G[1][v].end());
    }
    CHECK(V == vertex2id.size());
  }
  
  void CentralityNaive::PreCompute(const vector<pair<int, int> > &es, int) {
    BuildGraph(es);
    centrality_map = vector<double>(V, 0);
    const auto &forward_adj = G[0];
    
    for (int s = 0; (size_t)s < V; s++){
      // Brandes' algorithm (2001)
      // Step1. Compute distance and the number of shortest paths from s. O(m) time.
      vector<int>    distance(V, -1);
      vector<double> num_paths(V, 0);
      vector<double> delta(V, 0);
      vector<int>    order;
      queue<int>     que;
      
      distance[s] = 0;
      num_paths[s] = 1;
      que.push(s);
      
      while (!que.empty()){
        int v = que.front(); que.pop(); order.push_back(v);
        for (int w : forward_adj[v]){
          if (distance[w] == -1){
            distance[w] = distance[v] + 1;;
            que.push(w);
          }
          
          if (distance[w] == distance[v] + 1){
            num_paths[w] += num_paths[v];
          }
        }
      }
      
      // Step2. Aggregate betweenness centrality values. O(m) time.
      reverse(order.begin(), order.end());
      for (int v : order){
        for (int w : forward_adj[v]){
          if (distance[w] == distance[v] + 1 && num_paths[w] > 0){
            delta[v] += num_paths[v] / num_paths[w];
            delta[v] += num_paths[v] / num_paths[w] * delta[w];
          }
        }
      }
      for (int t = 0; (size_t)t < V; t++){
        if (s != t) centrality_map[t] += delta[t];
      }
    }
  }
}

