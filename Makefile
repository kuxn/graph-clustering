CC = g++
CXXFLAGS = -Wall -std=c++11 -O3

OBJECTS = main.o graph.o
TARGETS = main

all: $(TARGETS)

main: $(OBJECTS)
	$(CC) $(CXXFLAGS) -o $@ $(OBJECTS)

%.o:%.cpp
	$(CC) $(CXXFLAGS) -c $<

test: main
	./main

output: main
	time ./main > graph.dot

.PHONY: clean
clean:
	rm -f $(OBJECTS) $(TARGETS)
