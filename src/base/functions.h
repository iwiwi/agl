#include <string>
#include <sstream>


namespace agl {
// みさわさんに教えてもらったテク：
// template<typename ...> static inline int getchar_unlocked(void){ return getchar(); }

template <typename T, typename ...>
std::string to_string(const T& n) {
  std::ostringstream stm;
  stm << n;
  return stm.str();
}
}  // namespace agl
