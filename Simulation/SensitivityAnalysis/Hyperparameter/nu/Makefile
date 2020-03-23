CC=g++

all: BUSseq_nu1 BUSseq_nu2 BUSseq_nu3 BUSseq_nu4

# link
BUSseq_nu1: BUSseq_nu1.o rngstream.o
	$(CC) -fopenmp BUSseq_nu1.o rngstream.o -o BUSseq_nu1

BUSseq_nu2: BUSseq_nu2.o rngstream.o
	$(CC) -fopenmp BUSseq_nu2.o rngstream.o -o BUSseq_nu2

BUSseq_nu3: BUSseq_nu3.o rngstream.o
	$(CC) -fopenmp BUSseq_nu3.o rngstream.o -o BUSseq_nu3

BUSseq_nu4: BUSseq_nu4.o rngstream.o
	$(CC) -fopenmp BUSseq_nu4.o rngstream.o -o BUSseq_nu4



# compile
BUSseq_nu1.o: BUSseq_nu1.cpp
	$(CC) -c -fopenmp BUSseq_nu1.cpp -w -std=c++11

BUSseq_nu2.o: BUSseq_nu2.cpp
	$(CC) -c -fopenmp BUSseq_nu2.cpp -w -std=c++11

BUSseq_nu3.o: BUSseq_nu3.cpp
	$(CC) -c -fopenmp BUSseq_nu3.cpp -w -std=c++11

BUSseq_nu4.o: BUSseq_nu4.cpp
	$(CC) -c -fopenmp BUSseq_nu4.cpp -w -std=c++11

# clean all files
clean:
	rm -f *.o;
