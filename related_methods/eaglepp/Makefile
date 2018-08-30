CC=g++
CFLAGS=-c -Wall -I/usr/local/include/igraph -O3 -pipe  -I./include -ffast-math -ftree-vectorize 
LDFLAGS=-L/usr/local/lib -ligraph -fopenmp

CSOURCES=./src/Communities.cpp ./src/eagle.cpp ./src/MaxCache.cpp
ESOURCES=./src/main.cpp

EAGLE_OBJECTS=$(CSOURCES:.cpp=.o) $(ESOURCES:.cpp=.o)

EAGLE_NAME=bin/eagle

all: $(EAGLE_NAME)

$(EAGLE_NAME): $(EAGLE_OBJECTS)
	$(CC) $(LDFLAGS) $(EAGLE_OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm $(EAGLE_OBJECTS) $(EAGLE_NAME)

delbackups:
	rm *?.?*~

