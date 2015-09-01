#include "connectivity.h"
#include "gtest/gtest.h"
using namespace std;
using namespace agl;

TEST(is_connected, unweighted) {
  {
    unweighted_graph g({
      {0, 1},
      {1, 2},
    });
    ASSERT_TRUE(is_connected(g));
  }
  {
    unweighted_graph g({
      {0, 1},
      {2, 4},
    });
    ASSERT_FALSE(is_connected(g));
  }
}


TEST(is_connected, weighted) {
  {
    weighted_graph<double> g({
      {0, {1, 0.5}},
      {1, {2, 2.0}},
    });
    ASSERT_TRUE(is_connected(g));
  }
  {
    weighted_graph<double> g({
      {0, {1, 1.3}},
      {2, {4, 3.5}},
    });
    ASSERT_FALSE(is_connected(g));
  }
}
