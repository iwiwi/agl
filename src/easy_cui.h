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
#include "base/base.h"
#include "graph/graph.h"
using namespace std;
using namespace agl;

DEFINE_string(graph, "-", "input graph");

template<typename GraphType = G>
GraphType easy_cui_init(int argc, char **argv) {
  JLOG_INIT(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);
  GraphType g = read_graph_tsv<GraphType>(FLAGS_graph.c_str());
  pretty_print(g);
  return g;
}
