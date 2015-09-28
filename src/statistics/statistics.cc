#include "statistics.h"
#include "prettyprint.h"
using namespace std;

namespace agl {
pair<size_t, size_t> num_triangles_and_wedges(const G &g) {
  vector<bool> adj(g.num_vertices());

  size_t num_triangles = 0;
  size_t num_wedges = 0;
  for (V x : g.vertices()) {
    size_t dx = g.degree(x, kFwd) + g.degree(x, kBwd);

    for (V y : undirected_neighbors(g, x)) {
      size_t du = g.degree(y, kFwd) + g.degree(y, kBwd);
      if (tie(du, y) > tie(dx, x)) continue;

      for (V z : undirected_neighbors(g, y)) {
        if (adj[z]) ++num_triangles;
      }
      adj[y] = true;
    }

    size_t undirected_deg_x = 0;
    for (V u : undirected_neighbors(g, x)) {
      adj[u] = false;
      ++undirected_deg_x;
    }
    num_wedges += undirected_deg_x * (undirected_deg_x - 1) / 2;
  }

  num_wedges -= num_triangles * 3;
  return make_pair(num_triangles, num_wedges);
}

vector<pair<size_t, size_t>> num_local_triangles_and_wedges(const G &g) {
  vector<pair<size_t, size_t>> res(g.num_vertices());
  vector<bool> adj(g.num_vertices());

  for (V x : g.vertices()) {
    size_t dx = g.degree(x, kFwd) + g.degree(x, kBwd);

    for (V y : undirected_neighbors(g, x)) {
      size_t du = g.degree(y, kFwd) + g.degree(y, kBwd);
      if (tie(du, y) > tie(dx, x)) continue;

      for (V z : undirected_neighbors(g, y)) {
        if (adj[z]) {
          ++res[x].first;
          ++res[y].first;
          ++res[z].first;
        }
      }
      adj[y] = true;
    }

    size_t undirected_deg_x = 0;
    for (V u : undirected_neighbors(g, x)) {
      adj[u] = false;
      ++undirected_deg_x;
    }
    res[x].second = undirected_deg_x * (undirected_deg_x - 1) / 2;
  }

  for (V x : g.vertices()) {
    res[x].second -= res[x].first;
  }

  return res;
}

vector<double> local_clustering_coefficient(const G &g) {
  auto tw = num_local_triangles_and_wedges(g);
  vector<double> res(g.num_vertices());
  for (V v : g.vertices()) {
    if (tw[v].first + tw[v].second != 0) {
      res[v] = tw[v].first / (double)(tw[v].first + tw[v].second);
    }
  }
  return res;
}

double average_clustering_coefficient(const G &g) {
  auto lcc = local_clustering_coefficient(g);
  return accumulate(lcc.begin(), lcc.end(), 0.0) / g.num_vertices();
}

double global_clustering_coefficient(const G &g) {
  size_t num_triangles, num_wedges;
  tie(num_triangles, num_wedges) = num_triangles_and_wedges(g);
  size_t num_connected_triples = num_triangles + num_wedges;
  return num_triangles / (double)num_connected_triples;
}
}  // namespace agl
