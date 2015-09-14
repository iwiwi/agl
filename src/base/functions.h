#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include "type.h"

namespace agl {
// Cygwin/MinGW では to_string が定義されていないバグがある．
// http://stackoverflow.com/questions/12975341/to-string-is-not-a-member-of-std-says-so-g
// http://gcc.gnu.org/bugzilla/show_bug.cgi?id=52015
// これはみさわさんに教えてもらったテク：
// template<typename ...> static inline int getchar_unlocked(void){ return getchar(); }
template<typename T, typename ...>
std::string to_string(const T& n) {
  std::ostringstream stm;
  stm << n;
  return stm.str();
}

template<typename RangeType>
auto range_to_vector(RangeType range) -> std::vector<decltype(*range.begin())> {
  decltype(range_to_vector(range)) v;
  for (auto x : range) v.emplace_back(x);
  return v;
}

template<typename T>
std::vector<T> parse_space_separated_string(const std::string &str) {
  std::istringstream ss(str);
  std::vector<T> res;
  for (T t; ss >> t; ) res.emplace_back(t);
  return res;
}

template<typename T>
std::vector<T> parse_comma_separated_string(std::string str) {
  std::replace(str.begin(), str.end(), ',', ' ');
  return parse_space_separated_string<T>(str);
}

double get_current_time_sec();
}  // namespace agl
