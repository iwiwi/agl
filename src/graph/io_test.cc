#include "graph.h"
#include "gtest/gtest.h"
using namespace std;
using namespace agl;
using testing::Types;


namespace {
typedef Types<unweighted_graph, weighted_graph<int>, weighted_graph<double> > GraphTestTypes;

template<typename T>
class graph_io_test : public testing::Test {};
TYPED_TEST_CASE(graph_io_test, GraphTestTypes);

TYPED_TEST(graph_io_test, binary) {
  using GraphType = TypeParam;
  GraphType g(add_random_weight<GraphType>(generate_erdos_renyi(5, 2)));
  pretty_print(g);

  std::ostringstream oss (std::stringstream::binary);
  write_graph_binary(g, oss);

  std::istringstream iss(oss.str(), std::stringstream::binary);
  auto g2 = read_graph_binary<decltype(g)>(iss);

  CHECK(g.num_vertices() == g2.num_vertices());
  for(auto v : g.vertices()) {
    CHECK(g.degree(v) == g2.degree(v));
    for(size_t i = 0; i < g.degree(v); i++) {
      const auto& e1 = g.edge(v,i);
      const auto& e2 = g2.edge(v,i);
      CHECK(to(e1) == to(e2));
      CHECK(weight(e1) == weight(e2));
    }
  }
}

} // namespace


/*
namespace {

typedef Types<unweighted_graph, weighted_graph<int>> GraphTypes;
}

template<typename T>
class basic_graph_io_test : public testing::Test {};
TYPED_TEST_CASE(basic_graph_io_test, GraphTypes);

TYPED_TEST(basic_graph_io_test, tsv) {
  TypeParam g(add_weight<TypeParam>(gen_erdos_renyi(10, 2)));
  write_graph_tsv(g);
}
*/
