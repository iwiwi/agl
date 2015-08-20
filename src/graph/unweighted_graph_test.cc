#include "gtest/gtest.h"
#include "graph.h"
using namespace std;
using namespace agl;

TEST(unweighted_graph_test, instantiation) {
  G g;

  g.assign(unweighted_graph::edge_list_type{
          {0, {1}},
          {1, {2}},
          {2, {3}},
          {3, {1}},
  });

  pretty_print(g);
}

TEST(unweighted_graph_test, read_graph_tsv) {
  G g(gen_erdos_renyi(10, 2));
  pretty_print(g);
}

TEST(unweighted_graph_test, built_in) {
  G g = built_in_graph("karate_club");
  pretty_print(g);
}
