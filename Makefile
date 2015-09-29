.PHONY: all clean test

all:
	./waf configure && ./waf --targets="edge_centrality,edge_centrality_random_graph"

test:
	./waf configure && ./waf --targets="test" && bin/test

clean:
	./waf configure && ./waf clean
