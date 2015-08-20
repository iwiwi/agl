.PHONY: all clean

all:
	mkdir -p build && cd build && cmake .. && make -j

clean:
	rm -rf bin build
	