#include "graph.h"
#include "gtest/gtest.h"
using namespace std;
using namespace agl;
using testing::Types;

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
