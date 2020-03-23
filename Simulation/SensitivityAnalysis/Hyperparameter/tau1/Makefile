CC=g++

all: BUSseq_tau11 BUSseq_tau12 BUSseq_tau13 BUSseq_tau14

# link
BUSseq_tau11: BUSseq_tau11.o rngstream.o
	$(CC) -fopenmp BUSseq_tau11.o rngstream.o -o BUSseq_tau11

BUSseq_tau12: BUSseq_tau12.o rngstream.o
	$(CC) -fopenmp BUSseq_tau12.o rngstream.o -o BUSseq_tau12

BUSseq_tau13: BUSseq_tau13.o rngstream.o
	$(CC) -fopenmp BUSseq_tau13.o rngstream.o -o BUSseq_tau13

BUSseq_tau14: BUSseq_tau14.o rngstream.o
	$(CC) -fopenmp BUSseq_tau14.o rngstream.o -o BUSseq_tau14



# compile
BUSseq_tau11.o: BUSseq_tau11.cpp
	$(CC) -c -fopenmp BUSseq_tau11.cpp -w -std=c++11

BUSseq_tau12.o: BUSseq_tau12.cpp
	$(CC) -c -fopenmp BUSseq_tau12.cpp -w -std=c++11

BUSseq_tau13.o: BUSseq_tau13.cpp
	$(CC) -c -fopenmp BUSseq_tau13.cpp -w -std=c++11

BUSseq_tau14.o: BUSseq_tau14.cpp
	$(CC) -c -fopenmp BUSseq_tau14.cpp -w -std=c++11

# clean all files
clean:
	rm -f *.o;
