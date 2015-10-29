#include "vertex_centrality.h"
#include "shortest_path/shortest_path.h"
using namespace std;

DEFINE_int32(personalized_pagerank_num_iterations, 100, "");

namespace agl {
vector<double> vertex_centrality_degree(const G &g) {
  vector<double> vc(g.num_vertices());
  for (V v : g.vertices()) {
    vc[v] = g.degree(v);
  }
  return vc;
}

vector<double> vertex_centrality_classic_closeness(const G& g) {
  vector<double> cc(g.num_vertices());
  for (V v : g.vertices()) {
    auto ds = single_source_distance(g, v);
    for (V u : g.vertices()) cc[u] += ds[u];
  }
  for (auto &x : cc) x = 1 / x;
  return cc;
}

namespace {
// TODO: move to a proper position
vector<double> personalized_pagerank
(const G &g, const vector<pair<V, double>> &p, double C = 0.85) {
  const double kThreshold = false ? 0 : 1E-13;

  vector<double> x(g.num_vertices()), y(g.num_vertices());
  for (auto i : p) x[i.first] = i.second;

  for (int iter = 0; iter < FLAGS_personalized_pagerank_num_iterations; ++iter) {
    for (V u : g.vertices()) {
      if (g.degree(u) == 0) continue;

      double a = C * x[u] / g.degree(u);
      if (a > kThreshold) {
        for (V v : g.neighbors(u)) y[v] += a;
      }
      x[u] = 0;
    }
    for (auto i : p) x[i.first] += (1 - C) * i.second;
    x.swap(y);
  }
  return x;
}
}  // namespace

vector<double> vertex_centrality_pagerank(const G& g, double C) {
  vector<pair<V, double>> p(g.num_vertices(), {0, 1.0 / g.num_vertices()});
  for (V v : g.vertices()) p[v].first = v;
  return personalized_pagerank(g, p, C);
}
}  // namespace agl
