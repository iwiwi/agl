.PHONY: all test clean bin/test_all

all:
	mkdir -p build && cd build && cmake .. && make -j

bin/test_all:
	mkdir -p build && cd build && cmake .. && make test_all -j

test: bin/test_all
	bin/test_all

clean:
	rm -rf bin build
