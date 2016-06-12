CC 			= g++
CXXFLAGS 	= -Wall -std=c++11 -O3

OBJECTS 	= main.o graph.o lanczos.o
TARGET 		= main

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CXXFLAGS) -o $@ $(OBJECTS)

%.o:%.cpp
	$(CC) $(CXXFLAGS) -c $<

# Explicit dependencies required for headers
graph.o: 	graph.hpp
lanczos.o: 	lanczos.hpp


# Phony target to get around problem of having a file called 'clean'
.PHONY: clean
clean:
	rm -f $(OBJECTS) $(TARGET)

test: main
	./main

output: main
	time ./main > graph.dot

unit_tests: $(TARGET).o
	make -C unit_tests test

