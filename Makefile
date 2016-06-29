CXX 		= g++
CXXFLAGS	= -Wall -std=c++11 -O3 -finline-functions -ffast-math -fomit-frame-pointer -funroll-loops -DDebug
INCPATH		= include

OBJECTS 	= graph.o lanczos.o tqli.o partition.o analysis.o
TARGET 		= main
TESTTARGET 	= test

all: $(TARGET)

$(TARGET): $(OBJECTS) main.o
	$(CXX) $(CXXFLAGS) -I $(INCPATH) -o $@ $^

%.o:%.cc
	$(CXX) $(CXXFLAGS) -I $(INCPATH) -c $<

# Explicit dependencies required for headers
$(OBJECTS): 		$(INCPATH)/graph.h
test.o partition.o: $(INCPATH)/*.h lanczos.cc

define OBJECT_DEPENDS_ON_CORRESPONDING_HEADER
    $(1) : $(INCPATH)/${1:.o=.h}
endef

$(foreach object_file,$(OBJECTS),$(eval $(call OBJECT_DEPENDS_ON_CORRESPONDING_HEADER,$(object_file))))

# Phony target to get around problem of having a file called 'clean'
.PHONY: clean
clean:
	rm -f *.o $(TARGET) $(TESTTARGET)

$(TESTTARGET): $(OBJECTS) test.o
	$(CXX) $(CXXFLAGS) -I $(INCPATH) -o $@ $^
	time ./$(TESTTARGET)

output: main
	time ./main > graph.dot

eigenvalues_figure:
	gnuplot < eigenvalues.gnu

unit_tests: $(TARGET).o
	make -C unit_tests test
