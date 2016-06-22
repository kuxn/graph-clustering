CXX 		= g++
CXXFLAGS	= -Wall -std=c++11 -O3 -finline-functions -ffast-math -fomit-frame-pointer -funroll-loops
INCPATH		= include

all: t/wvtest

# EDIT HERE: add any "../*.c" files that you need to test
t/wvtest: wvtestmain.cc wvtest.cc t/wvtest.t.cc ../graph.cc ../lanczos.cc ../tqli.cc ../partition.cc ../test.cc ../analysis.cc
	$(CXX) $(CXXFLAGS) -D WVTEST_CONFIGURED -o $@ -I. -I../$(INCPATH) $^

runtests: all
	t/wvtest

test: all
	./wvtestrun $(MAKE) runtests

clean::
	rm -f *~ t/*~ *.o t/*.o t/wvtest
