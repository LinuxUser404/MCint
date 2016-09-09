


CC = g++
OPTS = -Wall -std=c++11 -O2 -DCONCURRENT=0
PROGS = mcint test_sob


INCLUDE = -I./include
SRC = -L.
LIBRARY =
MPI = -pthread -L/usr/lib64 -lmpi_cxx -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl

all: $(PROGS)

mcint:
	$(CC) $(OPTS) -o mcint mcint.cpp $(INCLUDE)
	
test_sob:
	$(CC) $(OPTS) -o test_sob test_sob.cpp $(INCLUDE) $(MPI)
	
clean:
	rm -f *.o $(PROGS)

test:
	./mcint -dim 12 -n 10000 -smooth 1 -seed 987654321
	./test_sob
