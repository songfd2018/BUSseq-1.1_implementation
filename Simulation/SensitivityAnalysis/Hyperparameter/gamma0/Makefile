CC=g++

all: BUSseq_gamma01 BUSseq_gamma02 BUSseq_gamma03 BUSseq_gamma04

# link
BUSseq_gamma01: BUSseq_gamma01.o rngstream.o
	$(CC) -fopenmp BUSseq_gamma01.o rngstream.o -o BUSseq_gamma01

BUSseq_gamma02: BUSseq_gamma02.o rngstream.o
	$(CC) -fopenmp BUSseq_gamma02.o rngstream.o -o BUSseq_gamma02

BUSseq_gamma03: BUSseq_gamma03.o rngstream.o
	$(CC) -fopenmp BUSseq_gamma03.o rngstream.o -o BUSseq_gamma03

BUSseq_gamma04: BUSseq_gamma04.o rngstream.o
	$(CC) -fopenmp BUSseq_gamma04.o rngstream.o -o BUSseq_gamma04



# compile
BUSseq_gamma01.o: BUSseq_gamma01.cpp
	$(CC) -c -fopenmp BUSseq_gamma01.cpp -w -std=c++11

BUSseq_gamma02.o: BUSseq_gamma02.cpp
	$(CC) -c -fopenmp BUSseq_gamma02.cpp -w -std=c++11

BUSseq_gamma03.o: BUSseq_gamma03.cpp
	$(CC) -c -fopenmp BUSseq_gamma03.cpp -w -std=c++11

BUSseq_gamma04.o: BUSseq_gamma04.cpp
	$(CC) -c -fopenmp BUSseq_gamma04.cpp -w -std=c++11

# clean all files
clean:
	rm -f *.o;
