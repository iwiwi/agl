#include "easy_cui.h"

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);
  printf("%.3f\t%.3f\n",
         average_clustering_coefficient(g),
         global_clustering_coefficient(g));
  return 0;
}
