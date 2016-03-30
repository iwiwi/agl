#include <easy_cui.h>

DEFINE_string(method, "gusfield", "gusfield, print_gomory_hu_tree");

#include "bi_dinitz.h"
#include "connected_components_filter.h"
#include "two_edge_CC_filter.h"
#include "cut_tree_with_2ECC.h"
#include "plain_gusfield.h"
#include "plain_gusfield_bi_dinitz.h"

G to_directed_graph(G&& g) {
  vector<pair<V, V>> ret;
  for (auto& e : g.edge_list()) {
    if (e.first < to(e.second)) ret.emplace_back(e.first, to(e.second));
    else if (to(e.second) < e.first) ret.emplace_back(to(e.second), e.first);
  }
  sort(ret.begin(), ret.end());
  ret.erase(unique(ret.begin(), ret.end()), ret.end());
  return G(ret);
}

using Gusfield = two_edge_CC_filter<cut_tree_with_2ECC>;
DEFINE_string(gomory_hu_builder, "cut_tree_with_2ECC", "cut_tree_with_2ECC");

string graph_name() {
  string x = FLAGS_graph;
  string ret;
  for (int i = int(x.size()) - 1; i >= 0; i--) {
    if (x[i] == '/' || x[i] == '\\') break;
    ret.push_back(x[i]);
  }
  reverse(ret.begin(), ret.end());
  return ret;
}

DEFINE_string(validation_data_path, "", "");

template<class gomory_hu_tree_t>
void print_gomory_hu_tree(G&& g) {

  fprintf(stderr, "print_gomory_hu_tree : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);

  if (FLAGS_validation_data_path == "") {
    auto gname = graph_name();
    FLAGS_validation_data_path = gname + ".tree";
  }
  gomory_hu_tree_t* gf;
  JLOG_PUT_BENCHMARK("gusfield_time") {
    gf = new gomory_hu_tree_t(g);
  }

  ofstream os(FLAGS_validation_data_path.c_str(), ios_base::out);
  gf->print_gomory_hu_tree(os);
  delete gf;
}

template<class T>
void main_(G&& g) {
  if (FLAGS_method == "gusfield") {
    JLOG_PUT_BENCHMARK("gusfield_time") {
      T gf(g);
    }
  } else if (FLAGS_method == "print_gomory_hu_tree") {
    print_gomory_hu_tree<T>(std::move(g));
  } else {
    fprintf(stderr, "unrecognized option '%s'\n", FLAGS_method.c_str());
    exit(-1);
  }

  JLOG_PUT("try_greedy_tree_packing", FLAGS_try_greedy_tree_packing);
  JLOG_PUT("getcap_counter", getcap_counter);
  JLOG_PUT("addcap_counter", addcap_counter);
  JLOG_PUT("preflow_eq_degree", preflow_eq_degree);
  JLOG_PUT("flow_eq_0", flow_eq_0);
  JLOG_PUT("gtp_edge_count_all", gtp_edge_count);
  JLOG_PUT("gtp_edge_miss_all", gtp_edge_miss);
  JLOG_PUT("gtp_edge_use_all", gtp_edge_use);
}

DEFINE_string(write_directed_graph_name,"", "");

int main(int argc, char** argv) {
  G g = easy_cui_init(argc, argv);
  fprintf(stderr, "easy_cui_init : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);
  if(FLAGS_graph.find(".directed") == string::npos) {
    g = to_directed_graph(std::move(g));
    fprintf(stderr, "load graph : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);
  }

  if(FLAGS_gomory_hu_builder == "write_directed_graph") {
    string output = FLAGS_write_directed_graph_name;
    if(output == "") {
      output = FLAGS_graph + ".directed"; 
    }
    write_graph_binary(g, output.c_str());
    exit(0);
  }

  if (FLAGS_gomory_hu_builder == "PlainGusfield") { 
    main_<plain_gusfield>(std::move(g));
  } else if (FLAGS_gomory_hu_builder == "PlainGusfield_bi_dinitz") { 
    main_<plain_gusfield_bi_dinitz>(std::move(g));
  } else if (FLAGS_gomory_hu_builder == "cut_tree_with_2ECC") {
    main_<Gusfield>(std::move(g));
  } else {
    fprintf(stderr, "unrecognized option -gomory_hu_builder='%s'\n", FLAGS_gomory_hu_builder.c_str());
    exit(-1);
  }
}
