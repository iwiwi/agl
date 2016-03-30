#pragma once
#include "connected_components_filter.h"

// 2ECC = 2-edge connected components
template<class handler_t>
class two_edge_CC_filter {
public:

  void lowlink_dfs(int v, int par, int& cur_ord) {
    lowlink[v] = order[v] = cur_ord++;
    for(int dir = 0; dir < 2; dir++) for (auto to : g.edges(v, D(dir))) {
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

  vector<vector<int>> get_local_id2global_id() const {
    const int cc_size = int(biconnected_graph_handler->handlers().size());
    vector<vector<int>> local_id2global_id(cc_size);

    const auto& local_indices = biconnected_graph_handler->local_indices();
    for(int v = 0; v < n_; v++) {
      auto& l2g = local_id2global_id[biconnected_graph_handler->handlers_index(v)];
      if (int(l2g.size()) < local_indices[v] + 1) l2g.resize(local_indices[v] + 1);
      l2g[local_indices[v]] = v;
    }

    return local_id2global_id;
  }

  two_edge_CC_filter(G& g) : n_(g.num_vertices()), g(g), uf(n_), lowlink(n_, -1), order(n_, -1) {

    G new_g;
    fprintf(stderr, "two_edge_CC_filter::constructor start : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);
    JLOG_ADD_BENCHMARK("time.decompose_2_connected_components") {
      for(int v = 0; v < n_; v++) for (auto& e : g.edges(v)) {
        uf.unite(v, to(e));
      }

      const int num_edges = g.num_edges();

      for(int v = 0; v < n_; v++) if (uf.root(v) == v) {
        int cur_ord = 0;
        lowlink_dfs(v, -1, cur_ord);
      }

      CHECK(bridge.size() + biconnected_graphs_edges.size() == size_t(num_edges));

      g.clear();
      new_g = G(biconnected_graphs_edges, n_);

      //dealloc
      lowlink.clear(); lowlink.shrink_to_fit();
      order.clear(); order.shrink_to_fit();
      // bridge.clear(); bridge.shrink_to_fit();
      biconnected_graphs_edges.clear(); biconnected_graphs_edges.shrink_to_fit();
    }
    fprintf(stderr, "two_edge_CC_filter::constructor end, before connected_components_filter : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);

    biconnected_graph_handler.reset(new connected_components_filter<handler_t>(new_g));
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

  void print_gomory_hu_tree(ostream& os) {
    vector<int> roots;
    for(int v = 0; v < n_; v++) if (uf.root(v) == v) roots.push_back(v);
    for(int i = 0; i < int(roots.size()) - 1; i++) os << roots[0] << " " << roots[i + 1] << " 0\n";
    for (auto& e : bridge) os << e.first << " " << e.second << " 1\n";

    vector<vector<int>> local_id2global_id = get_local_id2global_id();

    //weight2以上
    for(int i = 0; i < int(biconnected_graph_handler->handlers().size()); i++) {
      const auto& l2g = local_id2global_id[i];
      const auto& gusfield_core = biconnected_graph_handler->handlers()[i];
      for(int v = 0; v < int(gusfield_core.parent_weight().size()); v++) {
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
  const int n_;
  G& g;
  union_find uf;
  vector<int> lowlink, order;
  vector<pair<V, V>> bridge, biconnected_graphs_edges;

  unique_ptr<connected_components_filter<handler_t>> biconnected_graph_handler;
};