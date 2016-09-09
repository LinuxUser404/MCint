


CC = g++
OPTS = -Wall -std=c++11 -O2 -DCONCURRENT=0
PROGS = mcint test_sob


INCLUDE = -I./include
SRC = -L.
LIBRARY =

all: $(PROGS)

mcint:
	$(CC) $(OPTS) -o mcint mcint.cpp $(INCLUDE)
	
test_sob:

	
clean:
	rm -f *.o $(PROGS)

test:
	./mcint -dim 12 -n 10000 -smooth 1 -seed 987654321
#	./test_sob -n 0
