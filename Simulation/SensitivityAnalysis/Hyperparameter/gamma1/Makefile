CC=g++

all: BUSseq_gamma11 BUSseq_gamma12 BUSseq_gamma13 BUSseq_gamma14

# link
BUSseq_gamma11: BUSseq_gamma11.o rngstream.o
	$(CC) -fopenmp BUSseq_gamma11.o rngstream.o -o BUSseq_gamma11

BUSseq_gamma12: BUSseq_gamma12.o rngstream.o
	$(CC) -fopenmp BUSseq_gamma12.o rngstream.o -o BUSseq_gamma12

BUSseq_gamma13: BUSseq_gamma13.o rngstream.o
	$(CC) -fopenmp BUSseq_gamma13.o rngstream.o -o BUSseq_gamma13

BUSseq_gamma14: BUSseq_gamma14.o rngstream.o
	$(CC) -fopenmp BUSseq_gamma14.o rngstream.o -o BUSseq_gamma14



# compile
BUSseq_gamma11.o: BUSseq_gamma11.cpp
	$(CC) -c -fopenmp BUSseq_gamma11.cpp -w -std=c++11

BUSseq_gamma12.o: BUSseq_gamma12.cpp
	$(CC) -c -fopenmp BUSseq_gamma12.cpp -w -std=c++11

BUSseq_gamma13.o: BUSseq_gamma13.cpp
	$(CC) -c -fopenmp BUSseq_gamma13.cpp -w -std=c++11

BUSseq_gamma14.o: BUSseq_gamma14.cpp
	$(CC) -c -fopenmp BUSseq_gamma14.cpp -w -std=c++11

# clean all files
clean:
	rm -f *.o;
