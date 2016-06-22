CXX 		= g++
INCPATH		= include
CXXFLAGS 	= -Wall -std=c++11 -O3 -I $(INCPATH)

OBJECTS 	= main.o graph.o lanczos.o tqli.o partition.o
TESTOBJECTS = test.o graph.o lanczos.o tqli.o partition.o
TARGET 		= main
TESTTARGET 	= test

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS)

%.o:%.cc
	$(CXX) $(CXXFLAGS) -c $<

# Explicit dependencies required for headers
graph.o: 	$(INCPATH)/graph.h
lanczos.o: 	$(INCPATH)/graph.h $(INCPATH)/lanczos.h
tqli.o:		$(INCPATH)/tqli.h
patition.o:	$(INCPATH)/graph.h $(INCPATH)/lanczos.h $(INCPATH)/tqli.h $(INCPATH)/partition.h
test.o:		$(INCPATH)/graph.h $(INCPATH)/lanczos.h $(INCPATH)/tqli.h $(INCPATH)/partition.h $(INCPATH)/test.h
main.o:		$(INCPATH)/graph.h $(INCPATH)/partition.h

# Phony target to get around problem of having a file called 'clean'
.PHONY: clean
clean:
	rm -f $(OBJECTS) $(TARGET) $(TESTOBJECTS) $(TESTTARGET)

$(TESTTARGET): $(TESTOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(TESTOBJECTS)
	./$(TESTTARGET)

output: main
	time ./main > graph.dot

unit_tests: $(TARGET).o
	make -C unit_tests test

