#define _CRT_SECURE_NO_DEPRECATE

#include <fstream>
#include <cstdio>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <sstream>
#include <set>

using namespace std;
using ll = long long;

vector<char> fread(FILE* fp, int size) {
  vector<char> ret(size);
  fread(&ret[0], size, 1, fp);
  return ret;
}

template<typename T>
T fread(FILE* fp) {
  T ret;
  auto b = fread(&ret, sizeof(T), 1, fp);
  if(b != 1) {
    fprintf(stderr, "FAILED[*] ??.\n");
    exit(-1);
  }
  return ret;
}

void fwrite(FILE* fp, const string& text) {
  fwrite(&text[0], text.size(), 1, fp);
}

void fwrite(FILE* fp, int v) {
  fwrite(&v, sizeof(int), 1, fp);
}

void fwrite(FILE* fp, ll v) {
  fwrite(&v, sizeof(ll), 1, fp);
}

int main(int argc, char** argv) {

  if (argc < 3) {
    fprintf(stderr, "usage : graph.agl output.fam\n");
    exit(-1);
  }

  string input_file_name = argv[1];
  string output_file_name = argv[2];
  fprintf(stderr, "write fam to : '%s'.\n", output_file_name.c_str());
  FILE* ip = fopen(input_file_name.c_str(), "rb");

  const string kMagic = "AGL_BINARY";
  const string kVersion = "0.01";
  string header = kMagic + "\n" + kVersion + "\n" + "unweighted\n";
  auto h = fread(ip, header.size());
  string readheader(h.begin(), h.end());
  if (readheader != header) {
    fprintf(stderr, "FAILED[*] invalied header.\n");
  }

  int num_nodes = fread<int>(ip);
  ll num_edges = fread<ll>(ip);
  set<pair<int,int>> st;
  fprintf(stderr, "nodes = %d, edges = %lld\n", num_nodes, num_edges);
  for (int i = 0; i < num_nodes; i++) {
    ll deg = fread<ll>(ip);
    while (deg--) {
      int u = fread<int>(ip);
      if(i == u) continue;
      if(i < u) st.emplace(i, u);
      else st.emplace(u, i);
    }
  }

  fclose(ip);

  fprintf(stderr, "to undirected graph\n");
  fprintf(stderr, "nodes = %d, edges = %lld\n", num_nodes, st.size());
  FILE* op = fopen(output_file_name.c_str(), "wb");
  fprintf(op, "c agl to goldberg gam format.\n");
  fprintf(op, "p cut %d %lld\n", num_nodes, st.size());
  for(auto& uv : st) {
    fprintf(op, "a %d %d 1\n",uv.first + 1, uv.second + 1);
  }
  fclose(op);
  fprintf(stderr, "finished.\n");
}
