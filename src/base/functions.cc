#include "functions.h"

#include <sys/time.h>
#include <sys/utsname.h>
#include <unistd.h>
using namespace std;

namespace agl {
vector<string> split(const string &str, char splitter) {
  vector<string> res;
  string tmp;
  for (char c : str) {
    if (c == splitter) {
      if (tmp != "") res.emplace_back(move(tmp));
      tmp = "";
    } else {
      tmp += c;
    }
  }
  if (tmp != "") res.emplace_back(move(tmp));
  return res;
}

string strip(const string &s) {
  size_t i = 0;
  while (i < s.length() && isspace(s[i])) ++i;
  size_t j = s.length();
  while (i < j && isspace(s[j - 1])) --j;
  return s.substr(i, j - i);
}

double get_current_time_sec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}
}  // namespace agl
