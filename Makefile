CC = g++
CFLAGS = -Wall -std=c++11

OBJECTS = main.o graph.o
TARGETS = main

all: $(TARGETS)

main: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS)

%.o:%.cpp
	$(CC) $(CFLAGS) -c $<

test: main
	./main

.PHONY: clean
clean:
	rm -f $(OBJECTS) $(TARGETS)
