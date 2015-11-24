#include "box_cover.h"
using namespace std;

vector<V> box_cover_memb(const G &g, W radius) {
  V num_v = g.num_vertices();
  vector<vector<pair<V, W>>> node_lists;
  vector<pair<size_t, V>> center_candidates;
  for (V pv = 0; pv < num_v; ++pv) {
    queue<pair<V, W>> que;
    vector<bool> vis(num_v, false);
    vector<pair<V, W>> nodes;

    que.push(make_pair(pv, 0));
    nodes.push_back(make_pair(pv, 0));
    vis[pv] = true;

    while (!que.empty()) {
      V v = que.front().first;
      W dist = que.front().second;
      que.pop();
      if (dist >= radius) continue;
      for (V u : g.neighbors(v)) {
        if (vis[u]) continue;
        que.push(make_pair(u, dist + 1));
        nodes.push_back(make_pair(u, dist + 1));
        vis[u] = true;
      }
    }

    node_lists.push_back(nodes);
    center_candidates.push_back(make_pair(nodes.size(), pv));
  }

  map<size_t, set<V>> excluded_mass_of_non_centers;
  for (pair<size_t, V> p : center_candidates) {
    excluded_mass_of_non_centers[p.first].insert(p.second);
  }

  set<V> covered_nodes;
  set<V> center_nodes;
  vector<W> central_distance_of_node(num_v, num_v);
  while (covered_nodes.size() < (size_t)num_v) {
    V center_node_found;
    while (true) {
      V node, maximum_key;
      while (true) {
        map<size_t, set<V>>::iterator maximum_it =
            max_element(excluded_mass_of_non_centers.begin(),
                        excluded_mass_of_non_centers.end());
        maximum_key = maximum_it->first;
        set<V> &nodes = maximum_it->second;
        auto it = nodes.begin();
        advance(it, agl::random(nodes.size()));
        node = *it;
        if (center_nodes.find(node) != center_nodes.end()) {
          nodes.erase(it);
          if (nodes.empty()) excluded_mass_of_non_centers.erase(maximum_it);
        } else {
          break;
        }
      }
      map<V, W> nodes_visit;
      vector<pair<V, W>> &bfs = node_lists[node];
      V excluded_mass_of_node = 0;
      for (pair<V, W> p : bfs) {
        if (covered_nodes.find(p.first) == covered_nodes.end()) {
          excluded_mass_of_node++;
        }
      }
      excluded_mass_of_non_centers[maximum_key].erase(node);
      if (excluded_mass_of_non_centers[maximum_key].empty())
        excluded_mass_of_non_centers.erase(maximum_key);
      if (excluded_mass_of_node == maximum_key) {
        center_node_found = node;
        break;
      } else {
        excluded_mass_of_non_centers[excluded_mass_of_node].insert(node);
      }
    }
    center_nodes.insert(center_node_found);
    for (pair<V, W> p : node_lists[center_node_found]) {
      V i = p.first;
      W d = p.second;
      covered_nodes.insert(i);
      if (central_distance_of_node[i] > d) {
        central_distance_of_node[i] = d;
      }
    }
  }

  vector<int> ret(center_nodes.begin(), center_nodes.end());
  return ret;
}

vector<V> box_cover_burning(const G &g, W radius) { return {}; }