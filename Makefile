CXX 		= g++
#CXX 		= vtcxx -vt:c++ g++ -DVTRACE
CXXFLAGS	= -Wall -std=c++11 -O3 -finline-functions -ffast-math -fomit-frame-pointer -funroll-loops -g -DSO_
INCPATH		= include
VTINC		= /usr/local/include/vampirtrace
LIBDIR		= -L/usr/local/lib
LDLIBS		= -lboost_program_options

SRCDIR		= src
BUILDDIR	= build

CXXFILES	= $(shell find $(SRCDIR) -name '*.cc')		
OBJECTS 	= $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(CXXFILES:cc=o))
TARGET 		= main

all: $(TARGET)

$(TARGET): $(filter-out $(BUILDDIR)/test.o, $(OBJECTS))
	@echo "Linking..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBDIR) $(LDLIBS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cc
	@mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I $(INCPATH) -I $(VTINC) -c -o $@ $<

# Explicit dependencies required for headers
$(OBJECTS):	$(INCPATH)/graph.h
$(BUILDDIR)/test.o $(BUILDDIR)/partition.o:	$(INCPATH)/*.h $(SRCDIR)/lanczos.cc

define OBJECT_DEPENDS_ON_CORRESPONDING_HEADER
   $(1) : $(patsubst $(BUILDDIR)/%, $(INCPATH)/%, $(1:o=h)) 
endef

$(foreach object_file,$(OBJECTS),$(eval $(call OBJECT_DEPENDS_ON_CORRESPONDING_HEADER, $(filter-out $(BUILDDIR)/main.o, $(object_file)))))

# Phony target to get around problem of having a file called 'clean'
.PHONY: clean
clean:
	rm -rf *.z *.dSYM *.otf *.thumb $(BUILDDIR) $(TARGET)

tester: $(filter-out $(BUILDDIR)/main.o, $(OBJECTS))
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBDIR) $(LDLIBS)
	time ./tester

output: main
	time ./main > graph.dot

eigenvalues_figure:
	gnuplot < eigenvalues.gnu

unit_tests: $(BUILDDIR)/$(TARGET).o
	make -C unit_tests test
