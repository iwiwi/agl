#pragma once
#include "greedy_treepacking.h"

DECLARE_int32(try_greedy_tree_packing);
DECLARE_int32(try_large_degreepairs);
DECLARE_int32(separate_near_pairs_d);
DECLARE_int32(contraction_lower_bound);
DECLARE_bool(enable_greedy_tree_packing);
DECLARE_bool(enable_adjacent_cut);
DECLARE_bool(enable_goal_oriented_search);

namespace agl {

class disjoint_cut_set;
class separator;
class gomory_hu_tree_builder;

class cut_tree_with_2ECC {
  void find_cuts_by_tree_packing(std::vector<std::pair<V,V>>& edges, disjoint_cut_set* dcs, const std::vector<int>& degree);
  void contract_degree2_vertices(std::vector<std::pair<V,V>>& edges, std::vector<int>& degree);

  //次数の大きい頂点対をcutする
  void separate_high_degreepairs(separator* sep);

  //隣接頂点同士を見て、まだ切れていなかったらcutする
  void separate_adjacent_pairs(separator* sep);

  void separate_all(separator* sep);

  void separate_near_pairs(separator* sep);

  //次数の最も高い頂点に対して、出来る限りの頂点からflowを流してmincutを求める
  void find_cuts_by_goal_oriented_search(separator* sep);

public:

  cut_tree_with_2ECC(std::vector<std::pair<V, V>>&& edges, int num_vs);
  ~cut_tree_with_2ECC();

  int query(V u, V v) const;
  const std::vector<std::pair<V, int>>& parent_weight() const;

private:
  const int num_vertices_;
  std::unique_ptr<gomory_hu_tree_builder> gh_builder_;
};
} // namespace agl