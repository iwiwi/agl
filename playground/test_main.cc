#include <base/io.h>
#include <type_traits>
#include "easy_cui.h"

namespace {
	//set the end line code(regardless of the os)
	const string kEndLine = "\n";
	const string kMagic = "AGL_BINARY";
	const string kVersion = "0.01";

	// http://d.hatena.ne.jp/osyo-manga/20120211/1328922379
	extern void* enabler;
	template<bool B, typename T = void>
	using enabler_if = typename std::enable_if<B, T>::type*&;
}

namespace agl {

//format
template<typename GraphType> std::string graph_binary_format_core();
template<> string graph_binary_format_core<unweighted_graph>() {
	return "format=unweighted_graph";
}

template<typename GraphType>
const string graph_binary_format() {
	using decayed_gragh_type = typename std::decay<GraphType>::type;
	return graph_binary_format_core<decayed_gragh_type>();
}

// write_binary
template<typename EdgeType,
	enabler_if<std::is_pod<EdgeType>::value> = enabler >
void write_edge_binary(std::ostream& os, const EdgeType& e) {
	write_binary(os, e);
}

template<typename GraphType>
void write_graph_binary(const GraphType &g, std::ostream &os = std::cout) {
	//header
	os << kMagic << kEndLine << kVersion << kEndLine << graph_binary_format<GraphType>() << kEndLine;

	//body
	write_binary(os, g.num_vertices());
	write_binary(os, g.num_edges());
	for(V v = 0; v < g.num_vertices(); v++) {
		write_binary(os, g.degree(v));
		// todo : only unweighted edge
		for(const auto& e : g.neighbors(v)) {
				write_edge_binary(os, e);
		}
	}
  os.flush();
}


//reader binary
template<typename EdgeType,
	enabler_if<std::is_pod<EdgeType>::value> = enabler >
EdgeType read_edge_binary(std::istream& is) {
	unweighted_edge e;
	read_binary(is, &e);
	return e;
}

template<typename GraphType = G>
GraphType read_graph_binary(std::istream &is = std::cin) {
	//header
	std::string magic,version,format;
	std::getline(is, magic);
	std::getline(is, version);
	std::getline(is, format);

	CHECK_MSG(magic == kMagic, "Invalid file magic.");
	CHECK_MSG(version == kVersion, "Invalid file version.");
	CHECK_MSG(format == graph_binary_format<GraphType>(), "Invalid file format.");

	//body
	typename GraphType::V num_vertices;
	std::size_t num_edges;
	read_binary(is, &num_vertices);
	read_binary(is, &num_edges);

	std::vector<std::vector<typename GraphType::E>> edges(num_vertices);

	for(V v = 0; v < num_vertices; v++) {
		std::size_t degree;
		read_binary(is, &degree);
		edges[v].reserve(degree);
		for(std::size_t i = 0; i < degree; i++) {
			edges[v].emplace_back(read_edge_binary<typename GraphType::E>(is));
		}
	}

	GraphType deserialized_graph;
	deserialized_graph.assign(std::move(edges));
	return deserialized_graph;
}

} // namespace agl

DEFINE_string(output, "-", "output binary");

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);
  pretty_print(g);

  ofstream ofs(FLAGS_output, ios::out | ios::binary);
  write_graph_binary(g, ofs);
  ofs.close();

  ifstream ifs(FLAGS_output, ios::in | ios::binary);
  auto g2 = read_graph_binary<G>(ifs);

  for(auto v : g.vertices()) {
  	for(size_t i = 0; i < g.degree(v); i++) {
  		CHECK(g.edge(v,i) == g2.edge(v,i));
  	}
  }

  pretty_print(g2);
  return 0;
}
