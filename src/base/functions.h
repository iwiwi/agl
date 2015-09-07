#include <string>
#include <sstream>


namespace agl {
// Cygwin/MinGW では to_string が定義されていないバグがある．
// http://stackoverflow.com/questions/12975341/to-string-is-not-a-member-of-std-says-so-g
// http://gcc.gnu.org/bugzilla/show_bug.cgi?id=52015
// これはみさわさんに教えてもらったテク：
// template<typename ...> static inline int getchar_unlocked(void){ return getchar(); }

template <typename T, typename ...>
std::string to_string(const T& n) {
  std::ostringstream stm;
  stm << n;
  return stm.str();
}
}  // namespace agl
