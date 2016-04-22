#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>

#define CHECK(expr)                               \
  if (expr) {                                     \
  } else {                                        \
    fprintf(stderr, "CHECK Failed (%s:%d): %s\n", \
            __FILE__, __LINE__, #expr);           \
    exit(EXIT_FAILURE);                           \
  }

#define fst first
#define snd second

inline bool Equal(double a, double b){
  return abs(a - b) < 1e-9;
}

// for debug
template <typename S, typename T> std::ostream &operator<<(std::ostream &out, const std::pair<S, T> &p) {
  out << "(" << p.first << ", " << p.second << ")";
  return out;
}

template <typename T> std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  out << "[";
  for (size_t i = 0; i < v.size(); i++){
    if (i != 0) out << ", ";
    out << v[i];
  }
  out << "]";
  return out;
}

template<typename T> inline void SafeDelete(T* &p){
  if (p != nullptr){
    delete p;
    p = nullptr;
  }
}

template<typename T> inline void SafeDeleteArray(T* &p){
  if (p != nullptr){
    delete [] p;
    p = nullptr;
  }
}

#endif /* COMMON_H */
