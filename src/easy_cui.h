/*
 * 試作やベンチマーク目的でグラフを扱う CUI アプリケーションを作りたい場合，
 * こいつを使うと便利です．ズボラな人用．必ず本体から 1 度のみ include されるようにして下さい．
 */

#pragma once
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stack>
#include <queue>
#include <set>
#include <map>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include "agl.h"

DEFINE_string(type, "auto", "auto, tsv, agl, built_in, gen");
DEFINE_string(graph, "-", "input graph");
DEFINE_bool(force_undirected, false, "Automatically add reverse edges?");

using namespace std;
using namespace agl;

template<typename GraphType = G>
string guess_type() {
  FILE* fp = fopen(FLAGS_graph.c_str(), "r");
  if(fp == NULL) return "built_in";
  
  char buf[20];
  if(fgets(buf, sizeof(buf), fp) == NULL){
    string mes = "An error occured on read '" + FLAGS_graph + "'.";
    FAIL_MSG(mes.c_str());
  }
  fclose(fp);

  string header(buf);
  if(header == "AGL_BINARY\n") {
    return "agl";
  } else {
    return "tsv";
  }
}
