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
  fenwick_tree(size_t n) : x(n, 0) {}

  T sum(size_t l, size_t r) {
    if (l == 0) {
      T ans = 0;
      for (r; r >= 0; r = (r & (r + 1)) - 1) ans += x[r];
      return ans;
    }
    return sum(0, r) - sum(0, l - 1);
  }
  void add(size_t k, T a) {
    for (; k < x.size(); k |= k + 1) x[k] += a;
  }

 private:
  std::vector<T> x;
};
}  // namespace agl
