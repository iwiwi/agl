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
using namespace std;
using namespace agl;

DEFINE_string(type, "tsv", "tsv, built_in");
DEFINE_string(graph, "-", "input graph");
DEFINE_bool(force_undirected, false, "Automatically add reverse edges?");

template<typename GraphType = G>
GraphType easy_cui_init(int argc, char **argv) {
  JLOG_INIT(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);

  G::edge_list_type es;
  if (FLAGS_type == "tsv") {
    es = read_edge_list_tsv(FLAGS_graph.c_str());
  } else if (FLAGS_type == "built_in") {
    es = built_in_edge_list(FLAGS_graph.c_str());
  } else if (FLAGS_type == "gen") {
    istringstream iss(FLAGS_graph);
    string family;
    iss >> family;
    if (family == "barbell") {
      V n;
      if (!(iss >> n)) n = 4;
      es = generate_barbell(n);
    } else if (family == "grid") {
      size_t r, c;
      if (!(iss >> r)) r = 4;
      if (!(iss >> c)) c = r;
      es = generate_grid(r, c);
    } else if (family == "erdos_renyi") {
      V n;
      double d;
      if (!(iss >> n)) n = 10;
      if (!(iss >> d)) d = 3.0;
      es = generate_erdos_renyi(n, d);
    } else if (family == "random_planar") {
      V n;
      size_t e;
      if (!(iss >> n)) n = 10;
      if (!(iss >> e)) e = 25;
      es = generate_random_planar(n, e);
    } else if (family == "cycle") {
      V n;
      if (!(iss >> n)) n = 10;
      es = generate_cycle(n);
    } else if (family == "ba") {
      V n, m;
      if (!(iss >> n)) n = 10;
      if (!(iss >> m)) m = 5;
      es = generate_ba(n, m);
    } else if (family == "dms") {
      V n, m, k0;
      if (!(iss >> n)) n = 10;
      if (!(iss >> m)) m = 5;
      if (!(iss >> k0)) k0 = -2;
      es = generate_dms(n, m, k0);
    } else if (family == "hk") {
      V n, m;
      double p;
      if (!(iss >> n)) n = 10;
      if (!(iss >> m)) m = 5;
      if (!(iss >> p)) p = 0.5;
      es = generate_hk(n, m, p);
    } else if (family == "ws") {
      V n, d;
      double p;
      if (!(iss >> n)) n = 10;
      if (!(iss >> d)) d = 4;
      if (!(iss >> p)) p = 0.5;
      es = generate_ws(n, d, p);
    } else {
      FAIL_MSG("Unknown generator family: " + family);
    }
  }

  if (FLAGS_force_undirected) es = make_undirected(es);

  GraphType g(es);
  pretty_print(g);
  return g;
}
