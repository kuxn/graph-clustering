CXX 		= g++
INCPATH		= include
LIBDIR		= -L/usr/local/lib
LDLIBS		= -lboost_program_options
#GFLAGS		= -lprofiler
#CXXFLAGS	= -I $(INCPATH) -Wall -std=c++11 -O3 -finline-functions -ffast-math -fomit-frame-pointer -funroll-loops -DSO_ -DMedian_
CXXFLAGS	= -I $(INCPATH) -Wall -std=c++11 -O3 -finline-functions -ffast-math -fomit-frame-pointer -funroll-loops -DSO_ -DFiedler_

SRCDIR		= src
BUILDDIR	= build

CXXFILES	= $(shell find $(SRCDIR) -name '*.cc')		
OBJECTS 	= $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(CXXFILES:cc=o))
TARGET 		= ../main_serial

all: $(TARGET)

$(TARGET): $(filter-out $(BUILDDIR)/test.o, $(OBJECTS))
	@echo "Linking..."
	$(CXX) -o $@ $^ $(LIBDIR) $(LDLIBS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cc
	@mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) $(GFLAGS) -c -o $@ $<

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
	rm -rf $(BUILDDIR) $(TARGET)

serial_tester: $(filter-out $(BUILDDIR)/main.o, $(OBJECTS))
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBDIR) $(LDLIBS)
	mv serial_tester ../

unit_tests: $(BUILDDIR)/$(TARGET).o
	make -C unit_tests test
