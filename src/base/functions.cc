#include "functions.h"

#include <sys/time.h>
#include <sys/utsname.h>
#include <unistd.h>

namespace agl {
double get_current_time_sec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}
}  // namespace agl
