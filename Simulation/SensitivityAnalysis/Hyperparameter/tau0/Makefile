CC=g++

all: BUSseq_tau01 BUSseq_tau02 BUSseq_tau03 BUSseq_tau04

# link
BUSseq_tau01: BUSseq_tau01.o rngstream.o
	$(CC) -fopenmp BUSseq_tau01.o rngstream.o -o BUSseq_tau01

BUSseq_tau02: BUSseq_tau02.o rngstream.o
	$(CC) -fopenmp BUSseq_tau02.o rngstream.o -o BUSseq_tau02

BUSseq_tau03: BUSseq_tau03.o rngstream.o
	$(CC) -fopenmp BUSseq_tau03.o rngstream.o -o BUSseq_tau03

BUSseq_tau04: BUSseq_tau04.o rngstream.o
	$(CC) -fopenmp BUSseq_tau04.o rngstream.o -o BUSseq_tau04



# compile
BUSseq_tau01.o: BUSseq_tau01.cpp
	$(CC) -c -fopenmp BUSseq_tau01.cpp -w -std=c++11

BUSseq_tau02.o: BUSseq_tau02.cpp
	$(CC) -c -fopenmp BUSseq_tau02.cpp -w -std=c++11

BUSseq_tau03.o: BUSseq_tau03.cpp
	$(CC) -c -fopenmp BUSseq_tau03.cpp -w -std=c++11

BUSseq_tau04.o: BUSseq_tau04.cpp
	$(CC) -c -fopenmp BUSseq_tau04.cpp -w -std=c++11

# clean all files
clean:
	rm -f *.o;
