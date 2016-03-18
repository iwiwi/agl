#pragma once
#include <base/jlog.h>

namespace jlog_internal {

class jlog_conditional {
 public:
  jlog_conditional(bool add, const char *path, bool condition, bool glog = true)
      : path_(path), add_(add), condition_(condition), glog_(glog) {
    if(condition_) start_ = get_current_time_sec();
  }

  ~jlog_conditional() {
    if(condition_) {
    double r = get_current_time_sec() - start_;
      if (add_) {
        jlog::jlog_add(path_, r, glog_);
      } else {
        jlog::jlog_put(path_, r, glog_);
      }
    }
  }

  operator bool() {
    return false;
  }

 private:
  const char *path_;
  bool add_;
  bool condition_;
  double start_;
  bool glog_;
};

} // jlog_internal

#define JLOG_PUT_BENCHMARK_IF(...)                                  \
  if (jlog_internal::jlog_conditional o__                           \
      = jlog_internal::jlog_conditional(false, __VA_ARGS__)); else

#define JLOG_ADD_BENCHMARK_IF(...)                                  \
  if (jlog_internal::jlog_conditional o__                           \
      = jlog_internal::jlog_conditional(true, __VA_ARGS__)); else

