#include "graph_index_interface.h"
#include <unordered_set>
#include <memory>

namespace agl {
//
// Dynamic index update scenario.
// A scenario consists of
// (1) construction, and
// (2) update workloads (groups of update queries)
//
template<typename GraphType = G>
struct dynamic_index_evaluation_scenario {
  class update_interface;
  using edge_list_type = typename GraphType::edge_list_type;
  using workload_type = std::vector<std::shared_ptr<update_interface>>;
  using dynamic_index_type = graph_dynamic_index_interface<GraphType>;

  // Variables
  GraphType initial_graph;
  std::vector<workload_type> workloads;

  // Workload construction
  void add_workload_edge_addition(const edge_list_type &update_es);
  void add_workload_edge_addition_and_removal(const edge_list_type &update_es);
  void add_workload_edge_addition_and_removal_random(size_t num_update_es);

  // Evaluation
  struct workload_result;
  struct scenario_result;
  scenario_result evaluate(dynamic_index_type *idx);

  // Update queries
  class update_add_edge;
  class update_remove_edge;

 private:
  class nothing_index;
  std::vector<double> evaluate_internal(dynamic_index_type *idx);
};

//
// Definition of dynamic update queries
//
template<typename GraphType = G>
class dynamic_index_evaluation_scenario<GraphType>::update_interface {
 public:
  virtual void apply(GraphType *g, graph_dynamic_index_interface<GraphType> *i) const = 0;
};

template<typename GraphType = G>
class dynamic_index_evaluation_scenario<GraphType>::update_add_edge
: public dynamic_index_evaluation_scenario<GraphType>::update_interface {
 public:
  using E = typename GraphType::E;

  V v_from;
  E e;

  update_add_edge(V v_from, const E &e) : v_from(v_from), e(e) {}

  virtual void apply(GraphType *g, graph_dynamic_index_interface<GraphType> *i) const override {
    g->add_edge(v_from, e);
    i->add_edge(*g, v_from, e);
  }
};

template<typename GraphType = G>
class dynamic_index_evaluation_scenario<GraphType>::update_remove_edge
: public dynamic_index_evaluation_scenario<GraphType>::update_interface {
 public:
  V v_from, v_to;

  update_remove_edge(V v_from, V v_to) : v_from(v_from), v_to(v_to) {}

  virtual void apply(GraphType *g, graph_dynamic_index_interface<GraphType> *i) const override {
    g->remove_edge(v_from, v_to);
    i->remove_edge(*g, v_from, v_to);
  }
};

//
// Workload generation
//
template<typename GraphType>
void dynamic_index_evaluation_scenario<GraphType>::add_workload_edge_addition(const edge_list_type &update_es) {
  workloads.emplace_back();
  auto &w = workloads.back();
  w.resize(update_es.size());
  for (size_t i : make_irange(update_es.size())) {
    const auto &e = update_es[i];
    w[i] = std::make_shared<update_add_edge>(e.first, e.second);
  }
}

template<typename GraphType>
void dynamic_index_evaluation_scenario<GraphType>::add_workload_edge_addition_and_removal(const edge_list_type &update_es) {
  add_workload_edge_addition(update_es);

  workloads.emplace_back();
  auto &w = workloads.back();
  w.resize(update_es.size());
  for (size_t i : make_irange(update_es.size())) {
    const auto &e = update_es[update_es.size() - i - 1];
    w[i] = std::make_shared<update_remove_edge>(e.first, e.second);
  }
}

template<typename GraphType>
void dynamic_index_evaluation_scenario<GraphType>::add_workload_edge_addition_and_removal_random(size_t num_update_es) {
  edge_list_type update_es(num_update_es);
  std::unordered_set<std::pair<V, V>> update_es_set;

  std::uniform_int_distribution<V> rng(0, initial_graph.num_vertices() - 1);
  for (size_t i = 0; i < num_update_es; ++i) {
    V u, v;
    do {
      u = rng(agl::random);
      v = rng(agl::random);
    } while (u == v || is_adjacent(initial_graph, u, v) || update_es_set.count({u, v}));

    update_es_set.insert({u, v});
    update_es[i] = {u, v};
  }

  add_workload_edge_addition_and_removal(update_es);
}

//
// Workload evaluation
//
template<typename GraphType>
struct dynamic_index_evaluation_scenario<GraphType>::workload_result {
  double total_graph_and_index_udpate_time_sec;
  double total_graph_update_time_sec;  // Time used for just updating graph data structure
  double total_index_update_time_sec;  // Time used purely for index update
  double average_index_update_time_sec;
};

template<typename GraphType>
struct dynamic_index_evaluation_scenario<GraphType>::scenario_result {
  double construction_time_sec;
  std::vector<workload_result> workload_results;
};

template<typename GraphType>
typename dynamic_index_evaluation_scenario<GraphType>::scenario_result
dynamic_index_evaluation_scenario<GraphType>::evaluate(dynamic_index_type *i) {
  scenario_result r;
  {
    // Construction
    double t = get_current_time_sec();
    i->construct(initial_graph);
    r.construction_time_sec = get_current_time_sec() - t;
    JLOG_PUT("construction_time_sec", r.construction_time_sec);
  }
  {
    // Update workloads
    r.workload_results.resize(workloads.size());
    nothing_index ni;
    std::vector<double> t0 = evaluate_internal(&ni);
    std::vector<double> t1 = evaluate_internal(i);
    for (size_t i : make_irange(workloads.size())) {
      auto &wr = r.workload_results[i];
      wr.total_graph_and_index_udpate_time_sec = t1[i];
      wr.total_graph_update_time_sec = t0[i];
      wr.total_index_update_time_sec = t1[i] - t0[i];
      wr.average_index_update_time_sec = (t1[i] - t0[i]) / workloads[i].size();

      JLOG_ADD_OPEN("workloads") {
        JLOG_PUT("total_graph_and_index_udpate_time_sec", wr.total_graph_and_index_udpate_time_sec);
        JLOG_PUT("total_graph_update_time_sec", wr.total_graph_update_time_sec);
        JLOG_PUT("total_index_update_time_sec", wr.total_index_update_time_sec);
        JLOG_PUT("average_index_update_time_sec", wr.average_index_update_time_sec);
      }
    }
  }

  return r;
}

template<typename GraphType>
class dynamic_index_evaluation_scenario<GraphType>::nothing_index : public dynamic_index_type {
 public:
  virtual void construct(const GraphType &g) override {}
  virtual void add_edge(const GraphType &g, V v_from, const E &e) override {}
  virtual void remove_edge(const GraphType &g, V v_from, V v_to) override {}
  virtual void add_vertices(const GraphType &g, V old_num_vertices) override {}
  virtual void remove_vertices(const GraphType&, V) override {}
};

template<typename GraphType>
std::vector<double> dynamic_index_evaluation_scenario<GraphType>::evaluate_internal(dynamic_index_type *i) {
  std::vector<double> r;
  GraphType g = initial_graph;
  for (const auto &w : workloads) {
    double t = get_current_time_sec();
    for (const auto &u : w) {
      u->apply(&g, i);
    }
    r.emplace_back(get_current_time_sec() - t);
  }
  return r;
}

//
// Pretty print
//
template<typename GraphType>
void pretty_print(const dynamic_index_evaluation_scenario<GraphType> &s, std::ostream &os = std::cerr) {
  pretty_print(s.initial_graph, os);

  for (size_t i : make_irange(s.workloads.size())) {
    if (i > 0) os << "----------" << std::endl;
    os << "  Number of updates: " << s.workloads[i].size() << std::endl;
  }
  os << "=========" << std::endl;
}
}  // namespace agl
