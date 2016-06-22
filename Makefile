CXX 		= g++
CXXFLAGS	= -Wall -std=c++11 -O3 -finline-functions -ffast-math -fomit-frame-pointer -funroll-loops
INCPATH		= include

OBJECTS 	= main.o graph.o lanczos.o tqli.o partition.o analysis.o
TESTOBJECTS = test.o graph.o lanczos.o tqli.o partition.o analysis.o
TARGET 		= main
TESTTARGET 	= test

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -I $(INCPATH) -o $@ $(OBJECTS)

%.o:%.cc
	$(CXX) $(CXXFLAGS) -I $(INCPATH) -c $<

# Explicit dependencies required for headers
graph.o: 	$(INCPATH)/graph.h
lanczos.o: 	$(INCPATH)/graph.h $(INCPATH)/lanczos.h
tqli.o:		$(INCPATH)/tqli.h
patition.o:	$(INCPATH)/graph.h $(INCPATH)/lanczos.h $(INCPATH)/tqli.h $(INCPATH)/partition.h
test.o:		$(INCPATH)/graph.h $(INCPATH)/lanczos.h $(INCPATH)/tqli.h $(INCPATH)/partition.h $(INCPATH)/test.h
main.o:		$(INCPATH)/graph.h $(INCPATH)/partition.h
analysis.o:	$(INCPATH)/graph.h $(INCPATH)/partition.h

# Phony target to get around problem of having a file called 'clean'
.PHONY: clean
clean:
	rm -f $(OBJECTS) $(TARGET) $(TESTOBJECTS) $(TESTTARGET)

$(TESTTARGET): $(TESTOBJECTS)
	$(CXX) $(CXXFLAGS) -I $(INCPATH) -o $@ $(TESTOBJECTS)
	./$(TESTTARGET)

output: main
	time ./main > graph.dot

unit_tests: $(TARGET).o
	make -C unit_tests test

