#pragma once
#include "ConnectedComponentsFilter.h"

// 2ECC = 2-edge connected components
template<class handler_t>
class TwoEdgeCCFilter {
public:
  const int n;
  const G& g;

  void lowlink_dfs(int v, int par, int& cur_ord) {
    lowlink[v] = order[v] = cur_ord++;
    FOR(dir, 2) for (auto to : g.edges(v, D(dir))) {
      if (to == par) continue;
      if (order[to] == -1) {
        lowlink_dfs(to, v, cur_ord);
        lowlink[v] = min(lowlink[v], lowlink[to]);
        if (order[v] < lowlink[to]) bridge.emplace_back(v, to);
        else biconnected_graphs_edges.emplace_back(v, to);
      } else {
        lowlink[v] = min(lowlink[v], lowlink[to]);
        if (v < to) biconnected_graphs_edges.emplace_back(v, to);
      }
    }
  }

  TwoEdgeCCFilter(const G& g) : n(g.num_vertices()), g(g), uf(n), lowlink(n, -1), order(n, -1) {
    FOR(v, n) for (auto& e : g.edges(v)) {
      uf.unite(v, to(e));
    }

    const int num_edges = g.num_edges();

    FOR(v, n) if (uf.root(v) == v) {
      int cur_ord = 0;
      lowlink_dfs(v, -1, cur_ord);
    }

    CHECK(sz(bridge) + sz(biconnected_graphs_edges) == num_edges);

    G new_g(biconnected_graphs_edges, n);
    biconnected_graph_handler.reset(new ConnectedComponentsFilter<handler_t>(new_g));
  }

public:
  int query(V u, V v) {
    int ans = biconnected_graph_handler->query(u, v);
    if (ans == 0) {
      if (uf.is_same(u, v)) return 1; // 橋で間接的につながっている
      else return 0;
    }
    return ans;
  }

  void aggregate_gomory_hu_tree_weight() const {
    map<int,int> weight_num;
    const int wight0_num = biconnected_graph_handler->num_connected_components() - sz(bridge) - 1;
    weight_num[0] = wight0_num;
    const int wight1_num = sz(bridge);
    weight_num[1] = wight1_num;

    //weight2以上
    for(auto& gusfield_core : biconnected_graph_handler->handlers()) {
      for(auto& kv : gusfield_core.parent_weight()) {
        if(kv.first == -1) continue; // 親への辺が存在しない
        int weight = kv.second;
        CHECK(weight >= 2);
        weight_num[weight]++;
      }
    }

    for(const auto& wn : weight_num) {
      JLOG_ADD_OPEN("gomory-hu_edge") {
        JLOG_PUT("weight", wn.first);
        JLOG_PUT("count", wn.second);
      }
    }
  }

  void print_gomory_hu_tree(ostream& os) {
    const int n = sz(lowlink);
    vector<int> roots;
    FOR(v, n) if (uf.root(v) == v) roots.push_back(v);
    FOR(i, sz(roots) - 1) os << roots[0] << " " << roots[i + 1] << " 0\n";
    for (auto& e : bridge) os << e.first << " " << e.second << " 1\n";

    const int cc_size = sz(biconnected_graph_handler->handlers());
    vector<vector<int>> local_id2global_id(cc_size);

    const auto& local_indices = biconnected_graph_handler->local_indices();
    FOR(v, n) {
      auto& l2g = local_id2global_id[biconnected_graph_handler->handlers_index(v)];
      if (sz(l2g) < local_indices[v] + 1) l2g.resize(local_indices[v] + 1);
      l2g[local_indices[v]] = v;
    }

    //weight2以上
    FOR(i,sz(biconnected_graph_handler->handlers())) {
      const auto& l2g = local_id2global_id[i];
      const auto& gusfield_core = biconnected_graph_handler->handlers()[i];
      FOR(v, sz(gusfield_core.parent_weight())) {
        const auto& kv = gusfield_core.parent_weight()[v];
        int u = kv.first;
        if (u == -1) continue; // 親への辺が存在しない
        int weight = kv.second;
        CHECK(weight >= 2);
        os << l2g[v] << " " << l2g[u] << " " << weight << "\n";
      }
    }
  }

private:
  union_find uf;
  vector<int> lowlink, order;
  vector<pair<V, V>> bridge, biconnected_graphs_edges;

  unique_ptr<ConnectedComponentsFilter<handler_t>> biconnected_graph_handler;
};