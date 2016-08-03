CXX 		= g++
CXXFLAGS	= -Wall -std=c++11 -O3 -finline-functions -ffast-math -fomit-frame-pointer -funroll-loops -DSO_
INCPATH		= /serial/include
VTINC		= /usr/local/include/vampirtrace

all: t/wvtest
SRCDIR      = ../serial/src
CXXFILES  	= $(shell find $(SRCDIR) -name '*.cc')

# EDIT HERE: add any "../*.c" files that you need to test
t/wvtest: wvtestmain.cc wvtest.cc t/wvtest.t.cc $(filter-out $(SRCDIR)/main.cc, $(CXXFILES))
	$(CXX) $(CXXFLAGS) -D WVTEST_CONFIGURED -o $@ -I. -I../$(INCPATH) -I$(VTINC) $^

runtests: all
	t/wvtest

test: all
	./wvtestrun $(MAKE) runtests

clean::
	rm -f *~ t/*~ *.o t/*.o t/wvtest
