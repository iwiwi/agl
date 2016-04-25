#include "easy_cui.h"
#include "diameter.h"

int main(int argc, char** argv) {
  G g = easy_cui_init(argc, argv);
  
  JLOG_PUT_BENCHMARK("diameter") {
    printf("diameter = %d\n", diameter(g));
  }

  return 0;
}
