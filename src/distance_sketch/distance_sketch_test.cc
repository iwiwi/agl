#include "distance_sketch.h"
#include "gtest/gtest.h"
using namespace agl;
using namespace std;
using namespace agl::distance_sketch;

TEST(all_distances_sketches, hoge) {
  G g(gen_erdos_renyi(10, 2));
  auto rs = generate_rank_array(g.num_vertices());

  auto ads = compute_all_distances_sketches(g, 2, rs);
  auto srs = compute_sketch_retrieval_shortcuts(g, 2, rs);
  pretty_print(ads);
  pretty_print((const all_distances_sketches&)srs);

  for (V v : g.vertices()) {
    auto s = compute_all_distances_sketch_from(g, v, 2, rs);
    pretty_print(s);
    ASSERT_EQ(s, ads.retrieve_sketch(v));

    pretty_print(srs.retrieve_shortcuts(v));
    pretty_print(srs.retrieve_sketch(v));
    ASSERT_EQ(s, srs.retrieve_sketch(v));
  }
}
