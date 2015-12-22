/**
 * update_interface, update_add_edge, update_remove_edge = single update
 * update_workload = sequence of updates
 * update_scenario = initial graph + sequence of workloads
 */

#include <unordered_set>
#include <memory>
#include <algorithm>
#include "graph_index_interface.h"

namespace agl {
namespace dynamic_graph_update {
//
//! Single graph update
//
template<typename GraphType>
class update_interface {
 public:
  virtual ~update_interface() {}
  virtual void apply_to_graph(GraphType *g) const = 0;

  // The update should already have been applied to |g|.
  virtual void apply_to_index(GraphType *g, dynamic_graph_index_interface<GraphType> *i) const = 0;
};

template<typename GraphType>
class update_add_edge : public update_interface<GraphType> {
 public:
  using E = typename GraphType::E;
  V v_from;
  E e;

  update_add_edge(V v_from, const E &e) : v_from(v_from), e(e) {}
  virtual ~update_add_edge() {}

  virtual void apply_to_graph(GraphType *g) const override {
    g->add_edge(v_from, e);
  }

  virtual void apply_to_index(GraphType *g, dynamic_graph_index_interface<GraphType> *i) const override {
    i->add_edge(*g, v_from, e);
  }
};

template<typename GraphType>
class update_remove_edge : public update_interface<GraphType> {
 public:
  V v_from, v_to;

  update_remove_edge(V v_from, V v_to) : v_from(v_from), v_to(v_to) {}
  virtual ~update_remove_edge() {}

  virtual void apply_to_graph(GraphType *g) const override {
    g->remove_edge(v_from, v_to);
  }

  virtual void apply_to_index(GraphType *g, dynamic_graph_index_interface<GraphType> *i) const override {
    i->remove_edge(*g, v_from, v_to);
  }
};

//
//! Update workload (a sequence of graph update)
//
template<typename GraphType>
using update_workload = std::vector<std::shared_ptr<update_interface<GraphType>>>;

//! Measure the total time consumption for updating the graph and an index
template<typename GraphType>
double evaluate_workload(GraphType *g, dynamic_graph_index_interface<GraphType> *i, const update_workload<GraphType> &w) {
  assert(g != nullptr);
  double t = get_current_time_sec();
  for (const auto &u : w) {
    u->apply_to_graph(g);
    if (i != nullptr) u->apply_to_index(g, i);
  }
  return get_current_time_sec() - t;
}

//
//! Update scenario (initial graph + graph update workloads)
//
template<typename GraphType>
struct update_scenario {
  using edge_list_type = typename GraphType::edge_list_type;

  V initial_num_vertices;
  edge_list_type initial_edges;
  std::vector<update_workload<GraphType>> update_workloads;
};

//
// Evaluation using |update_scenario|
//

//! A struct to represent the evaluation result
struct scenario_evaluation_result {
  struct workload_evaluation_result {
    double total_graph_and_index_udpate_time_sec;
    double total_graph_update_time_sec;  // Time used for just updating graph data structure
    double total_index_update_time_sec;  // Time used purely for index update
    double average_index_update_time_sec;
  };

  double construction_time_sec;
  std::vector<workload_evaluation_result> workload_results;
};

//! Evaluate a dynamic indexing method by using the whole scenario
template<typename GraphType>
scenario_evaluation_result evaluate_scenario
(dynamic_graph_index_interface<GraphType> *i, const update_scenario<GraphType> &s) {
  scenario_evaluation_result r;
  std::vector<double> t0, t1;
  {
    GraphType g1(s.initial_edges);
    {
      // Construction
      double t = get_current_time_sec();
      i->construct(g1);
      r.construction_time_sec = get_current_time_sec() - t;
      JLOG_PUT("construction_time_sec", r.construction_time_sec);
    }
    for (const auto &w : s.update_workloads) {
      // Update
      t1.emplace_back(evaluate_workload(&g1, i, w));
    }
  }
  {
    GraphType g0(s.initial_edges);
    for (const auto &w : s.update_workloads) {
      t0.emplace_back(evaluate_workload(&g0, (dynamic_graph_index_interface<GraphType>*)nullptr, w));
    }
  }

  std::cout << s.update_workloads.size() << std::endl;

  // Aggregation
  r.workload_results.resize(s.update_workloads.size());
  for (size_t i : make_irange(s.update_workloads.size())) {
    auto &wr = r.workload_results[i];
    wr.total_graph_and_index_udpate_time_sec = t1[i];
    wr.total_graph_update_time_sec = t0[i];
    wr.total_index_update_time_sec = t1[i] - t0[i];
    wr.average_index_update_time_sec = (t1[i] - t0[i]) / s.update_workloads[i].size();

    JLOG_ADD_OPEN("workloads") {
      JLOG_PUT("total_graph_and_index_udpate_time_sec", wr.total_graph_and_index_udpate_time_sec);
      JLOG_PUT("total_graph_update_time_sec", wr.total_graph_update_time_sec);
      JLOG_PUT("total_index_update_time_sec", wr.total_index_update_time_sec);
      JLOG_PUT("average_index_update_time_sec", wr.average_index_update_time_sec);
    }
  }

  return r;
}

//
// Scenario generation
//
template<typename GraphType = G>
update_scenario<GraphType> generate_scenario_random_addition_and_removal
(const typename GraphType::edge_list_type &edges, size_t num_update_es) {
  update_scenario<GraphType> us;
  us.initial_edges.assign(edges.begin(), edges.end());
  auto &es = us.initial_edges;
  std::sort(es.begin(), es.end());
  V num_vs = us.initial_num_vertices = num_vertices_from_edge_list(es);

  us.update_workloads.assign(2, update_workload<GraphType>(num_update_es));
  std::unordered_set<std::pair<V, V>> update_es_set;
  std::uniform_int_distribution<V> rng(0, num_vs - 1);
  for (size_t i = 0; i < num_update_es; ++i) {
    V u, v;
    for (;;) {
      u = rng(agl::random);
      v = rng(agl::random);
      if (u == v) continue;
      if (std::binary_search(es.begin(), es.end(), std::make_pair(u, v))) continue;
      if (update_es_set.count({u, v})) continue;
      break;
    }

    update_es_set.insert({u, v});
    us.update_workloads[0][i] =
        std::make_shared<update_add_edge<GraphType>>(u, v);
    us.update_workloads[1][num_update_es - i - 1] =
        std::make_shared<update_remove_edge<GraphType>>(u, v);
  }

  return us;
}
}  // namespace dynamic_graph_update
}  // namespace agl
