// Copyright 2013, Takuya Akiba
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Takuya Akiba nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "jlog.h"

std::string FLAGS_jlog_out = "./jlog";
bool FLAGS_jlog_suppress_log = false;
std::string FLAGS_jlog_file_name = "";

namespace jlog_internal {
jlog jlog::instance_;

struct null_streambuf : public std::streambuf {
  virtual int overflow(int c) { return c; }
};

null_streambuf ns;
std::ostream jlog::null_ostream(&ns);

double get_current_time_sec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

long get_memory_usage() {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return usage.ru_maxrss;
}
}

void JLOG_INIT(int *argc, char **argv) {
  int k = 0;
  for (int i = 0; i < *argc; ++i) {
    if (strncmp(argv[i], "--jlog_out=", 11) == 0) {
      FLAGS_jlog_out = argv[i] + 11;
    } else if (strncmp(argv[i], "--jlog_suppress_log", 19) == 0) {
      FLAGS_jlog_suppress_log = true;
    } else if(strncmp(argv[i], "--jlog_file_name=", 17) == 0) {
      FLAGS_jlog_file_name = argv[i] + 17;
    } else {
      argv[k++] = argv[i];
    }
  }
  *argc = k;
  jlog_internal::jlog::init(*argc, argv);
}

