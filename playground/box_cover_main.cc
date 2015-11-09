#include "easy_cui.h"
#include "box_cover.h"

namespace {
// What portion of the vertices are really covered by the set?
double coverage(const G &g, vector<V> s, W rad) {
  // TODO: implement me
  return 0.0;
}
}

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);
  CHECK_MSG(FLAGS_force_undirected, "undirected only!!!");

  vector<pair<string, function<vector<V>(const G&, W)>>> algos{
    {"MEMB", box_cover_memb},
    {"Schneider", box_cover_burning},
  };

  for (auto a : algos) {
    JLOG_ADD_OPEN("algorithms") {
      JLOG_PUT("name", a.first);
      auto f = a.second;

      for (W rad = 0; rad <= 10; ++rad) {
        vector<V> res;
        JLOG_ADD_BENCHMARK("time") res = f(g, rad);
        JLOG_ADD("size", res.size());
        JLOG_PUT("coverage", coverage(g, res, rad));
      }
    }
  }
}
