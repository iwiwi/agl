#include "distance_sketch.h"
#include "gtest/gtest.h"
#include <iostream>
using namespace agl;
using namespace std;
using namespace agl::distance_sketch;

TEST(distance_sketch, hoge) {
  for (int trial = 0; trial < 100; ++trial) {
    G g(generate_erdos_renyi(10, 2));
    auto rs = generate_rank_array(g.num_vertices());

    auto ads = compute_all_distances_sketches(g, 2, rs);
    auto srs1 = compute_sketch_retrieval_shortcuts_via_ads_naive(g, 2, rs);
    auto srs2 = compute_sketch_retrieval_shortcuts_via_ads_fast(g, 2, rs);
    auto srs3 = compute_sketch_retrieval_shortcuts_via_ads_unweighted(g, 2, rs);
    pretty_print(ads);
    pretty_print((const all_distances_sketches&)srs1);

    for (V v : g.vertices()) {
      auto s = compute_all_distances_sketch_from(g, v, 2, rs);
      pretty_print(s);
      ASSERT_EQ(s, ads.retrieve_sketch(g, v));

      pretty_print(srs1.retrieve_shortcuts(g, v));
      pretty_print(srs1.retrieve_sketch(g, v));
      ASSERT_EQ(s, srs1.retrieve_sketch(g, v));
      ASSERT_EQ(s, srs2.retrieve_sketch(g, v));
      ASSERT_EQ(s, srs3.retrieve_sketch(g, v));
    }
  }
}

TEST(distance_sketch, shortcuts) {
  for (int trial = 0; trial < 100; ++trial) {
    G g(generate_erdos_renyi(10, 2));
    auto rs = generate_rank_array(g.num_vertices());
    // graphviz_draw_graph(g, "tmp.png");

    auto srs1 = compute_sketch_retrieval_shortcuts_via_ads_naive(g, 2, rs);
    auto srs2 = compute_sketch_retrieval_shortcuts_via_ads_fast(g, 2, rs);
    auto srs3 = compute_sketch_retrieval_shortcuts_via_ads_unweighted(g, 2, rs);

    for (V v : g.vertices()) {
      cout << "--" << endl;
      pretty_print(srs1.retrieve_shortcuts(g, v));
      pretty_print(srs1.retrieve_sketch(g, v));
      pretty_print(srs2.retrieve_shortcuts(g, v));
      pretty_print(srs2.retrieve_sketch(g, v));
      ASSERT_EQ(srs1.retrieve_shortcuts(g, v), srs2.retrieve_shortcuts(g, v));
      ASSERT_EQ(srs1.retrieve_shortcuts(g, v), srs3.retrieve_shortcuts(g, v));
    }
  }
}

TEST(distance_sketch, update) {
  static constexpr V kNumVertices = 20;
  static constexpr size_t kK = 2;

  for (int trial = 0; trial < 10; ++trial) {
    dynamic_index_evaluation_scenario<G> s;
    s.initial_graph.assign(generate_erdos_renyi(kNumVertices, 2));
    s.add_workload_edge_addition_and_removal_random(30);

    auto rs = generate_rank_array(kNumVertices);

    G g = s.initial_graph;
    dynamic_all_distances_sketches dads(kK, rs);
    dads.construct(g);

    for (const auto &w : s.workloads) {
      for (const auto &u : w) {
        //graphviz_draw_graph(g, "before.png");
        pretty_print(g);
        u->apply(&g, &dads);
        //graphviz_draw_graph(g, "after.png");
        pretty_print(g);
        auto ads = compute_all_distances_sketches(g, kK, rs);
        for (V v : g.vertices()) {
          cout << v << endl;
          pretty_print(ads.retrieve_sketch(g, v));
          pretty_print(dads.retrieve_sketch(g, v));
          ASSERT_EQ(dads.retrieve_sketch(g, v), ads.retrieve_sketch(g, v));
        }
      }
      pretty_print(g);
    }
  }
}