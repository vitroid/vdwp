CXX=g++
CXXFLAGS=-std=c++14 -O -Iopt-3.19/src #-g#-static #-fopenmp
LDFLAGS=opt-3.19/src/libopt.a #-lopt
BIN=histogram

%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)
all: $(BIN)

# get opt library
opt-3.19:
	wget https://www.decompile.com/download/not_invented_here/opt-3.19.tar.gz
	tar zxvf opt-3.19.tar.gz
opt-3.19/src/libopt.a: opt-3.19
	cd opt-3.19 && ./configure && make

clean:
	-rm histogram
