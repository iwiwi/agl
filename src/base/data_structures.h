#pragma once

#include <cstdint>
#include <type_traits>
#include <vector>

namespace agl {
class union_find {
 public:
  using T = int32_t;

  explicit union_find(T n = 0) {
    init(n);
  }

  void clear() {
    init(0);
  }

  void init(T n) {
    par_.assign(n, ~0);
  }

  T root(T x) {
    return par_[x] < 0 ? x : par_[x] = root(par_[x]);
  }

  void unite(T x, T y) {
    x = root(x);
    y = root(y);
    if (x == y) return;

    if (rank(x) < rank(y)) {
      par_[x] = y;
    } else {
      if (rank(x) == rank(y)) {
        par_[x] = ~(rank(x) + 1);
      }
      par_[y] = x;
    }
  }

  bool is_same(T x, T y) {
    return root(x) == root(y);
  }

 private:
  static_assert(std::is_signed<T>::value, "T should be signed");
  static_assert(std::is_integral<T>::value, "T should be integer");

  std::vector<T> par_;

  T rank(T x) {
    return ~par_[x];
  }
};

template<class T>
class fenwick_tree {
public:
  fenwick_tree(std::size_t n) : x(n, 0) {}

  T sum(std::size_t l, std::size_t r) {
    if (l == 0) {
      T ans = 0;
      for (int64_t k = r; k >= 0; k = (k & (k + 1)) - 1) ans += x[k];
      return ans;
    }
    return sum(0, r) - sum(0, l - 1);
  }

  void add(std::size_t k, T a) {
    for (; k < x.size(); k |= k + 1) x[k] += a;
  }

  /** Returns the index i of the first element such that sum(0, i) is
   *  not less than val. Returns the size of the tree if no such i exists.
   *  Note that sum(0, i) must be non-decreasing, i.e., all the elements
   *  in the tree must be non-negative.
   *  Time complexity: O(log n).
   */
  std::size_t lower_bound(T val) {
    int64_t n = x.size();
    if (n == 0) return 0;
    int64_t index = 1;
    int64_t width = 1;

    while (index <= n) {
      index *= 2;
      width *= 2;
    }
    index--;

    std::size_t ans = n;
    T sum = 0;

    while (width) {
      width /= 2;
      if (index >= n) {
        index -= width;
      } else {
        sum += x[index];
        if (sum >= val) {
          ans = index;
          sum -= x[index];
          index -= width;
        } else {
          index += width;
        }
      }
    }

    return ans;
  }

  /** Returns the index i of the first element such that sum(0, i) is
   *  greater than val. Returns the size of the tree if no such i exists.
   *  Note that sum(0, i) must be non-decreasing, i.e., all the elements
   *  in the tree must be non-negative.
   *  Time complexity: O(log n).
   */
  std::size_t upper_bound(T val) {
    int64_t n = x.size();
    if (n == 0) return 0;
    int64_t index = 1;
    int64_t width = 1;

    while (index <= n) {
      index *= 2;
      width *= 2;
    }
    index--;

    std::size_t ans = n;
    T sum = 0;

    while (width) {
      width /= 2;
      if (index >= n) {
        index -= width;
      } else {
        sum += x[index];
        if (sum > val) {
          ans = index;
          sum -= x[index];
          index -= width;
        } else {
          index += width;
        }
      }
    }

    return ans;
  }

 private:
  std::vector<T> x;
};
}  // namespace agl
