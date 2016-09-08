


CC = mpic++
OPTS = -Wall -O2
PROGS = mcint test_sob


INCLUDE = -I./include
SRC = -L.
LIBRARY =

all: $(PROGS)

mcint:
	$(CC) $(OPTS)    -o mcint mcint.cc $(INCLUDE) $(SRC) $(LIBRARY)
	
test_sob:
	$(CC) $(OPTS)    -o test_sob test_sob.cc $(INCLUDE) $(SRC) $(LIBRARY)
	
clean:
	rm -f *.o $(PROGS)

test:
	./mcint -dim 12 -n 10000 -smooth 1 -seed 987654321

