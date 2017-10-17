CC = gcc
CFLAGS = -W -Wall
LDFLAGS = -fopenmp
executables = sieve.out gauss.out

$(executable): $(objects)
	$(CC) $< $(LIB) -o $@

.PHONY: all clean fullclean sieve-test gauss-test graph

all: sieve.out gauss.out

clean:
	rm -f $(executables)

fullclean: clean
	rm -f *.dat *.ps

sieve.out: sieve.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

gauss.out: gauss.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

sieve-test: sieve.out
	./sieve.out

gauss-test: gauss.out
	./gauss.out

graph: sieve.out gauss.out plot.gnu
	./sieve.out -p
	./gauss.out -p
	gnuplot plot.gnu
