#include "data_structures.h"
#include "base/random.h"
#include "gtest/gtest.h"
#include <vector>
#include <algorithm>

using namespace std;
using namespace agl;

TEST(fenwick_tree_test, random_case) {
  for (int trial = 0; trial < 10; ++trial) {
    int n = 1 + agl::random(1000);
    fenwick_tree<int> t(n);
    int naive[n];
    
    for (int i = 0; i < n; ++i) {
      int a = agl::random(1000);
      t.add(i, a);
      naive[i] = a;
    }

    for (int i = 0; i < n; ++i) {
      ASSERT_EQ(t.sum(0, i), accumulate(naive, naive + i + 1, 0));
    }

    int sum = t.sum(0, n - 1);

    for (int i = 0; i < 1000; ++i) {
      int a = agl::random(sum) + 1;
      int index = t.lower_bound(a);
      ASSERT_TRUE(index == 0 || accumulate(naive, naive + index, 0) < a);
      ASSERT_TRUE(index == n || accumulate(naive, naive + index + 1, 0) >= a);
    }

    for (int i = 0; i < 1000; ++i) {
      int a = agl::random(sum) + 1;
      int index = t.upper_bound(a);
      ASSERT_TRUE(index == 0 || accumulate(naive, naive + index, 0) <= a);
      ASSERT_TRUE(index == n || accumulate(naive, naive + index + 1, 0) > a);
    }
  }
}
