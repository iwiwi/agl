#include <easy_cui.h>

DEFINE_string(s1, "s1.txt", "gomory-hu tree 1");
DEFINE_string(s2, "s2.txt", "gomory-hu tree 2");

#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())

map<int, vector<pair<V, V>>> load(const string& path) {
  FILE* fp = fopen(path.c_str(), "r");
  map<int, vector<pair<V, V>>> ret;
  int s, v, weight;
  while (~fscanf(fp, "%d%d%d", &s, &v, &weight)) {
    ret[weight].emplace_back(s, v);
  }
  return ret;
}

void check(map<int, vector<pair<V, V>>>& l, map<int, vector<pair<V, V>>>& r) {
  using ull = unsigned long long;
  int n = 1;
  for (auto& kv : l) n += sz(kv.second);
  vector<ull> zobrist_hash(n);
  FOR(i, n) zobrist_hash[i] = (ull(agl::random()) << 32) | ull(agl::random());

  union_find ufl(n), ufr(n);
  vector<ull> hl = zobrist_hash, hr = zobrist_hash;

  vector<int> weight;
  for (auto& kv : l) weight.push_back(kv.first);
  reverse(weight.begin(), weight.end());

  auto on_mismatched = [&](V u) {
    set<int> onlyl,onlyr,intersect;
    FOR(a, n) if (ufl.is_same(u, a)) onlyl.insert(a);
    FOR(a, n) if (ufr.is_same(u, a)) {
      if(onlyl.count(a)) onlyl.erase(a), intersect.insert(a);
      else onlyr.insert(a);
    } 
    int cnt;
    printf("left: ");
    cnt = 0;
    for(auto x : onlyl) {
      printf("%d, ",x);
      cnt++;
      if(cnt >= 10) break;
    }
    puts("");
    printf("right: ");
    cnt = 0;
    for(auto x : onlyr) {
      printf("%d, ",x);
      cnt++;
      if(cnt >= 10) break;
    }
    puts("");
    printf("intersect: ");
    cnt = 0;
    for(auto x : intersect) {
      printf("%d, ",x);
      cnt++;
      if(cnt >= 10) break;
    }
    puts("");
  };

  for (const int w : weight) {
    if(sz(l[w]) != sz(r[w])) printf("*");
    printf("w = %d, L : %d , R : %d", w, sz(l[w]), sz(r[w]));
    puts("");
  }

  for (const int w : weight) {
    // CHECK(sz(l[w]) == sz(r[w]));
    for (auto& uv : l[w]) {
      int u, v; tie(u, v) = uv;
      u = ufl.root(u), v = ufl.root(v);
      CHECK(u != v);
      ufl.unite(u, v);
      ull new_hash = hl[u] ^ hl[v];
      hl[ufl.root(u)] = new_hash;
    }

    for (auto& uv : r[w]) {
      int u, v; tie(u, v) = uv;
      u = ufr.root(u), v = ufr.root(v);
      CHECK(u != v);
      ufr.unite(u, v);
      ull new_hash = hr[u] ^ hr[v];
      hr[ufr.root(u)] = new_hash;
    }

    for (auto& uv : l[w]) {
      int u, v; tie(u, v) = uv;
      if (hl[ufl.root(u)] != hr[ufr.root(u)]) {
        on_mismatched(u);
        exit(-1);
      }
      if (hl[ufl.root(v)] != hr[ufr.root(v)]) {
        on_mismatched(v);
        exit(-1);
      }
    }
    for (auto& uv : r[w]) {
      int u, v; tie(u, v) = uv;
      if (hl[ufl.root(u)] != hr[ufr.root(u)]) {
        on_mismatched(u);
        exit(-1);
      }
      if (hl[ufl.root(v)] != hr[ufr.root(v)]) {
        on_mismatched(v);
        exit(-1);
      }
    }
  }
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  auto gomory_hu1 = load(FLAGS_s1);
  auto gomory_hu2 = load(FLAGS_s2);
  check(gomory_hu1, gomory_hu2);

  fprintf(stderr, "ok. %s == %s\n", FLAGS_s1.c_str(), FLAGS_s2.c_str());
  return 0;
}
