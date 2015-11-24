#include "box_cover.h"
using namespace std;

int main() {
  map<size_t, set<V>> excluded_mass_of_non_centers;
  set<V> center_nodes;
  excluded_mass_of_non_centers[10].insert(3);
  excluded_mass_of_non_centers[10].insert(5);
  excluded_mass_of_non_centers[20].insert(7);
  excluded_mass_of_non_centers[20].insert(9);
  excluded_mass_of_non_centers[30].insert(2);

  map<size_t, set<V>>::iterator maximum_it = max_element(
      excluded_mass_of_non_centers.begin(), excluded_mass_of_non_centers.end());
  set<V> &nodes = maximum_it->second;
  auto it = nodes.begin();
  advance(it, agl::random(nodes.size()));
  V node = *it;
  center_nodes.insert(node);
  if (center_nodes.find(node) != center_nodes.end()) {
    cerr << "it" << node << endl;
    nodes.erase(it);
    if (nodes.empty()) {
      excluded_mass_of_non_centers.erase(maximum_it);
    }
  }

  for (auto a : excluded_mass_of_non_centers) {
    cerr << a.first << endl;
    for (auto b : a.second) {
      cerr << b << endl;
    }
    cerr << endl;
  }
}