CXX 		= mpic++
#CXX 		= vtcxx -vt:c++ mpic++ -DVTRACE
CXXFLAGS	= -Wall -std=c++11 -O3 -finline-functions -ffast-math -fomit-frame-pointer -funroll-loops -DSO_
VTINC		= /usr/local/include/vampirtrace
LIBDIR		= -L/usr/local/lib
LDLIBS		= -lboost_mpi -lboost_serialization -lboost_program_options

INCPATH		= ../parallel/include
SRCDIR		= ../parallel/src
BUILDDIR	= ../parallel/build

CXXFILES	= $(shell find $(SRCDIR) -name '*.cc')		
OBJECTS 	= $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(CXXFILES:cc=o))
TARGET 		= main_parallel

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
	rm -rf data times output sbatch* slurm*

plot:
	chmod +x ./plot.sh
	./plot.sh

bench: $(TARGET)
	chmod +x ./submit.sh
	./submit.sh
