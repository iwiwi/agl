#include "distance_sketch.h"
#include "gtest/gtest.h"
#include <iostream>
using namespace agl;
using namespace std;
using namespace agl::distance_sketch;
using namespace agl::dynamic_graph_update;

TEST(distance_sketch, entry_comparison) {
  ASSERT_GT(entry(1, 2).raw, entry(2, 1).raw);
}

TEST(distance_sketch, ads_static_small) {
  for (int trial = 0; trial < 3; ++trial) {
    G g(generate_erdos_renyi(10, 2));
    // graphviz_draw_graph(g, "tmp.png");
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
      pretty_print(ads.retrieve_sketch(g, v));
      ASSERT_EQ(s, ads.retrieve_sketch(g, v));

      pretty_print(srs1.retrieve_shortcuts(g, v));
      pretty_print(srs1.retrieve_sketch(g, v));
      ASSERT_EQ(s, srs1.retrieve_sketch(g, v));
      ASSERT_EQ(s, srs2.retrieve_sketch(g, v));
      ASSERT_EQ(s, srs3.retrieve_sketch(g, v));
    }
  }
}


TEST(distance_sketch, ads_static_large) {
  for (int k : {1, 4, 16}) {
    for (int trial = 0; trial < 100; ++trial) {
      G g(generate_erdos_renyi(100, 2));
      auto rs = generate_rank_array(g.num_vertices());
      auto ads = compute_all_distances_sketches(g, k, rs);
      auto srs1 = compute_sketch_retrieval_shortcuts_via_ads_naive(g, k, rs);
      auto srs2 = compute_sketch_retrieval_shortcuts_via_ads_fast(g, k, rs);
      auto srs3 = compute_sketch_retrieval_shortcuts_via_ads_unweighted(g, k, rs);

      for (V v : g.vertices()) {
        auto s = compute_all_distances_sketch_from(g, v, k, rs);
        ASSERT_EQ(s, ads.retrieve_sketch(g, v));
        ASSERT_EQ(s, srs1.retrieve_sketch(g, v));
        ASSERT_EQ(s, srs2.retrieve_sketch(g, v));
        ASSERT_EQ(s, srs3.retrieve_sketch(g, v));
      }
    }
  }
}

TEST(distance_sketch, srs_static_small) {
  for (int trial = 0; trial < 3; ++trial) {
    for (int k : {2 /*1, 4, 16*/}) {
      G g(generate_erdos_renyi(100, 2));
      auto rs = generate_rank_array(g.num_vertices());
      // graphviz_draw_graph(g, "tmp.png");

      auto srs1 = compute_sketch_retrieval_shortcuts_via_ads_naive(g, k, rs);
      auto srs2 = compute_sketch_retrieval_shortcuts_via_ads_fast(g, k, rs);
      auto srs3 = compute_sketch_retrieval_shortcuts_via_ads_unweighted(g, k, rs);

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
}

TEST(distance_sketch, srs_static_large) {
  for (int k : {1, 4, 16}) {
    for (int trial = 0; trial < 100; ++trial) {
      G g(generate_erdos_renyi(100, 2));
      auto rs = generate_rank_array(g.num_vertices());
      auto srs1 = compute_sketch_retrieval_shortcuts_via_ads_naive(g, k, rs);
      auto srs2 = compute_sketch_retrieval_shortcuts_via_ads_fast(g, k, rs);
      auto srs3 = compute_sketch_retrieval_shortcuts_via_ads_unweighted(g, k, rs);

      for (V v : g.vertices()) {
        ASSERT_EQ(srs1.retrieve_shortcuts(g, v), srs2.retrieve_shortcuts(g, v));
        ASSERT_EQ(srs1.retrieve_shortcuts(g, v), srs3.retrieve_shortcuts(g, v));
      }
    }
  }
}


TEST(distance_sketch, ads_update_small) {
  static constexpr V kNumVertices = 10;
  static constexpr size_t kK = 2;

  for (int trial = 0; trial < 10; ++trial) {
    auto es = generate_erdos_renyi(kNumVertices, 2);
    auto us = generate_scenario_random_addition_and_removal(es, 10);
    us.initial_num_vertices = kNumVertices;

    auto rs = generate_rank_array(kNumVertices);

    G g(us.initial_edges, kNumVertices);
    dynamic_all_distances_sketches dads(kK, rs);
    dads.construct(g);

    for (const auto &w : us.update_workloads) {
      for (const auto &u : w) {
        //graphviz_draw_graph(g, "before.png");
        pretty_print(g);
        u->apply_to_graph(&g);
        u->apply_to_index(&g, &dads);
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

TEST(distance_sketch, ads_update_large) {
  static constexpr V kNumVertices = 100;

  for (int k : {1, 4, 16}) {
    for (int trial = 0; trial < 100; ++trial) {
      auto es = generate_erdos_renyi(kNumVertices, 2);
      auto us = generate_scenario_random_addition_and_removal(es, 100);
      us.initial_num_vertices = kNumVertices;

      auto rs = generate_rank_array(kNumVertices);

      G g(us.initial_edges, kNumVertices);
      dynamic_all_distances_sketches dads(k, rs);
      dads.construct(g);

      for (const auto &w : us.update_workloads) {
        for (const auto &u : w) {
          u->apply_to_graph(&g);
          u->apply_to_index(&g, &dads);
        }

        auto ads = compute_all_distances_sketches(g, k, rs);
        for (V v : g.vertices()) {
          ASSERT_EQ(dads.retrieve_sketch(g, v), ads.retrieve_sketch(g, v));
        }
      }
    }
  }
}


#if 0 ///////////////////////////////////////////////////////////////////////
TEST(distance_sketch, srs_update_small) {
  static constexpr V kNumVertices = 10;

  for (int k : {2 /*1, 4, 16*/}) {
    for (int trial = 0; trial < 1000; ++trial) {
      cout << "TRIAL: " << trial << endl;

      auto es = generate_random_planar(kNumVertices, kNumVertices);
      auto us = generate_scenario_random_addition_and_removal(es, 10);
      us.initial_num_vertices = kNumVertices;

      auto rs = generate_rank_array(kNumVertices);

      G g1(us.initial_edges, kNumVertices);
      G g2(us.initial_edges, kNumVertices);
      dynamic_all_distances_sketches dads(k, rs);
      dynamic_sketch_retrieval_shortcuts dsrs(k, rs);
      dads.construct(g1);
      dsrs.construct(g2);

      for (const auto &w : us.update_workloads) {
        for (const auto &u : w) {
          u->apply_to_graph(&g1);
          u->apply_to_index(&g1, &dads);
        }

        //  puts("nya");
       //   pretty_print(g2);
       //   pretty_print(compute_all_distances_sketch_from(g2, 1, k, rs));
      //    printf("1: "); pretty_print(dsrs.retrieve_sketch(g2, 1));
      //    printf("3: "); pretty_print(dsrs.retrieve_sketch(g2, 3));
        for (const auto &u : w) {
          //graphviz_draw_graph(g2, "0.png");
          u->apply_to_graph(&g2);
          u->apply_to_index(&g2, &dsrs);
          //graphviz_draw_graph(g2, "1.png");

          auto ads = compute_all_distances_sketches(g2, k, rs);
          auto srs = compute_sketch_retrieval_shortcuts(g2, k, rs);

          pretty_print(g2);
          pretty_print(compute_all_distances_sketches(g2, k, rs));
          pretty_print((const all_distances_sketches&)compute_sketch_retrieval_shortcuts(g2, k, rs));
          pretty_print((const all_distances_sketches&)dsrs.srs_);

          for (V v : g1.vertices()) {
            ASSERT_EQ(dsrs.retrieve_shortcuts(g2, v), srs.retrieve_shortcuts(g2, v)) << v;
            pretty_print(dsrs.retrieve_sketch(g2, v));
            pretty_print(ads.retrieve_sketch(g2, v));
            ASSERT_EQ(dsrs.retrieve_sketch(g2, v), ads.retrieve_sketch(g2, v)) << v;
          }

        }

        auto ads = compute_all_distances_sketches(g1, k, rs);
        auto srs = compute_sketch_retrieval_shortcuts(g2, k, rs);
        for (V v : g1.vertices()) {
          //pretty_print(g1);
          //cout << v << endl;
       //   pretty_print(ads.retrieve_sketch(g1, v));
       //   pretty_print(dads.retrieve_sketch(g1, v));
       //   pretty_print(dsrs.retrieve_sketch(g1, v));
          ASSERT_EQ(dads.retrieve_sketch(g1, v), ads.retrieve_sketch(g1, v));
          // ASSERT_EQ(dsrs.retrieve_shortcuts(g2, v), srs.retrieve_shortcuts(g2, v));
          ASSERT_EQ(dsrs.retrieve_sketch(g1, v), ads.retrieve_sketch(g1, v));
        }

        break;
      }
    }
  }
}

TEST(distance_sketch, srs_update_large) {
  static constexpr V kNumVertices = 100;

  for (int trial = 0; trial < 1000; ++trial) {
    for (int k : {1, 4, 16}) {
      auto es = generate_erdos_renyi(kNumVertices, 2);
      auto us = generate_scenario_random_addition_and_removal(es, 10);
      us.initial_num_vertices = kNumVertices;

      auto rs = generate_rank_array(kNumVertices);

      G g(us.initial_edges, kNumVertices);
      dynamic_sketch_retrieval_shortcuts dsrs(k, rs);
      dsrs.construct(g);

      for (const auto &w : us.update_workloads) {
        for (const auto &u : w) {
          u->apply_to_graph(&g);
          u->apply_to_index(&g, &dsrs);
        }

        auto ads = compute_all_distances_sketches(g, k, rs);
        auto srs = compute_sketch_retrieval_shortcuts(g, k, rs);
        for (V v : g.vertices()) {
          ASSERT_EQ(dsrs.retrieve_shortcuts(g, v), srs.retrieve_shortcuts(g, v));
          ASSERT_EQ(dsrs.retrieve_sketch(g, v), ads.retrieve_sketch(g, v));
        }

        break;
      }
    }
  }
}
#endif ///////////////////////////////////////////////////////////////////////