#pragma once
#include "agl.h"
#include "distance_sketch.h"

DECLARE_int32(distance_sketch_k);

namespace agl {
namespace distance_sketch {
///////////////////////////////////////////////////////////////////////////////
// Sketch Entry
///////////////////////////////////////////////////////////////////////////////
static const int kNumBitsVertex = 24;
static const int kNumBitsDistance = 8;

union entry {
  struct {
    int v : kNumBitsVertex;
    int d : kNumBitsDistance;
  };
  uint32_t raw;

  entry(V v = 0, int d = 0) {
    this->v = v;
    this->d = d;
  }

  bool operator==(const entry &e) const {
    return raw == e.raw;
  }

  bool operator<(const entry &e) const {
    return std::make_pair(d, v) < std::make_pair(e.d, e.v);
  }

  static constexpr V max_vertices() { return (1 << kNumBitsVertex) - 1; }
  static constexpr int max_distance() { return (1 << kNumBitsDistance) - 1; }
};

inline void pretty_print(const entry &e, std::ostream &ofs = std::cerr) {
  ofs << "(" << e.v << ", " << e.d << ")";
}

// TODO: weighted entry type


///////////////////////////////////////////////////////////////////////////////
// Rank
///////////////////////////////////////////////////////////////////////////////
using rank_type = uint64_t;
using rank_array = std::vector<rank_type>;

rank_array generate_rank_array(V num_vertices);

void pretty_print_rank_array(const rank_array &rank_array, std::ostream &ofs = std::cerr);

///////////////////////////////////////////////////////////////////////////////
// Static Sketch
///////////////////////////////////////////////////////////////////////////////
using vertex_sketch_raw = std::vector<entry>;

class graph_sketches_interface {
 public:
  virtual ~graph_sketches_interface() {}

  // Access
  virtual vertex_sketch_raw retrieve_sketch(const G &g, V v) = 0;

  // Statistics
  virtual double average_size() const = 0;
};

struct all_distances_sketches : public graph_sketches_interface {
  size_t k;
  std::vector<vertex_sketch_raw> sketches;
  rank_array ranks;

  virtual vertex_sketch_raw retrieve_sketch(const G &g, V v) override {
    return sketches[v];
  }

  virtual double average_size() const override {
    auto f = [](double r, const vertex_sketch_raw &s) -> double {
      return r + s.size();
    };
    return std::accumulate(sketches.begin(), sketches.end(), 0.0, f) / sketches.size();
  }

  virtual ~all_distances_sketches() {}
};

// Remove unnecessary entries
vertex_sketch_raw purify_sketch(vertex_sketch_raw sketch, size_t k, const rank_array &ranks);

vertex_sketch_raw compute_all_distances_sketch_from
(const G &g, V v, size_t k, const rank_array &ranks, D d = kFwd);

all_distances_sketches compute_all_distances_sketches
(const G &g, size_t k = FLAGS_distance_sketch_k, const rank_array &ranks = {}, D d = kFwd);

void pretty_print(const vertex_sketch_raw &s, std::ostream &ofs = std::cerr);
void pretty_print(const all_distances_sketches &ads, std::ostream &ofs = std::cerr);

///////////////////////////////////////////////////////////////////////////////
// Sketch Retrieval Shortcuts
///////////////////////////////////////////////////////////////////////////////
struct sketch_retrieval_shortcuts : public all_distances_sketches {
  // TODO: neighbor removal
  virtual vertex_sketch_raw retrieve_sketch(const G &g, V v) override;

  virtual vertex_sketch_raw retrieve_shortcuts(const G &g, V v) {
    return all_distances_sketches::retrieve_sketch(g, v);
  }
};

sketch_retrieval_shortcuts compute_sketch_retrieval_shortcuts
(const G &g, size_t k = FLAGS_distance_sketch_k, const rank_array &ranks = {}, D d = kFwd);

sketch_retrieval_shortcuts compute_sketch_retrieval_shortcuts_via_ads_naive
(const G &g, size_t k = FLAGS_distance_sketch_k, const rank_array &ranks = {}, D d = kFwd);

sketch_retrieval_shortcuts compute_sketch_retrieval_shortcuts_via_ads_fast
(const G &g, size_t k = FLAGS_distance_sketch_k, const rank_array &ranks = {}, D d = kFwd);

sketch_retrieval_shortcuts compute_sketch_retrieval_shortcuts_via_ads_unweighted
(const G &g, size_t k = FLAGS_distance_sketch_k, const rank_array &ranks = {}, D d = kFwd);

///////////////////////////////////////////////////////////////////////////////
// Dynamization
//////////////////////////////////////////////////////////////////////////////
class dynamic_graph_sketches : public dynamic_graph_index_interface<G> {
 public:
  virtual vertex_sketch_raw retrieve_sketch(const G &g, V v) = 0;
  virtual double average_sketch_length() const = 0;
};

class dynamic_all_distances_sketches : public dynamic_graph_sketches {
 public:
  dynamic_all_distances_sketches
  (size_t k = FLAGS_distance_sketch_k, const rank_array &ranks = {}, D d = kFwd) :
    d_(d) {
    ads_.k = k;
    ads_.ranks = ranks;
  }

  virtual void construct(const G &g) override {
    ads_ = compute_all_distances_sketches(g, ads_.k, ads_.ranks, d_);
  }

  virtual vertex_sketch_raw retrieve_sketch(const G &g, V v) override {
    return ads_.retrieve_sketch(g, v);
  }

  virtual double average_sketch_length() const override {
    return ads_.average_size();
  }

  virtual void add_edge(const G &g, V v_from, const E &e) override;

  virtual void remove_edge(const G &g, V v_from, V v_to) override;

  virtual void add_vertices(const G &g, V old_num_vertices) override {
    assert(false);
  }

  virtual void remove_vertices(const G&, V) override {
    assert(false);
  }

 private:
  D d_;
  all_distances_sketches ads_;

  // Helpers for update
  bool add_entry(V v, V s, W d);
};

class dynamic_sketch_retrieval_shortcuts : public dynamic_graph_sketches {
 public:
  dynamic_sketch_retrieval_shortcuts
  (size_t k = FLAGS_distance_sketch_k, const rank_array &ranks = {}, D d = kFwd) :
    d_(d) {
    srs_.k = k;
    srs_.ranks = ranks;
  }

  virtual void construct(const G &g) override {
    srs_ = compute_sketch_retrieval_shortcuts(g, srs_.k, srs_.ranks, d_);
  }

  virtual vertex_sketch_raw retrieve_sketch(const G &g, V v) override {
    return srs_.retrieve_sketch(g, v);
  }

  virtual double average_sketch_length() const override {
    return srs_.average_size();
  }

  virtual void add_edge(const G &g, V v_from, const E &e) override {}

  virtual void remove_edge(const G &g, V v_from, V v_to) override {}

  virtual void add_vertices(const G &g, V old_num_vertices) override {
    assert(false);
  }

  virtual void remove_vertices(const G&, V) override {
    assert(false);
  }

 private:
  D d_;
  sketch_retrieval_shortcuts srs_;
};
}  // namespace distance_sketch
}  // namespace agl
