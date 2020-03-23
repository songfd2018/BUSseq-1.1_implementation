CC=g++

all: BUSseq_alpha1 BUSseq_alpha2 BUSseq_alpha3 BUSseq_alpha4

# link
BUSseq_alpha1: BUSseq_alpha1.o rngstream.o
	$(CC) -fopenmp BUSseq_alpha1.o rngstream.o -o BUSseq_alpha1

BUSseq_alpha2: BUSseq_alpha2.o rngstream.o
	$(CC) -fopenmp BUSseq_alpha2.o rngstream.o -o BUSseq_alpha2

BUSseq_alpha3: BUSseq_alpha3.o rngstream.o
	$(CC) -fopenmp BUSseq_alpha3.o rngstream.o -o BUSseq_alpha3

BUSseq_alpha4: BUSseq_alpha4.o rngstream.o
	$(CC) -fopenmp BUSseq_alpha4.o rngstream.o -o BUSseq_alpha4



# compile
BUSseq_alpha1.o: BUSseq_alpha1.cpp
	$(CC) -c -fopenmp BUSseq_alpha1.cpp -w -std=c++11

BUSseq_alpha2.o: BUSseq_alpha2.cpp
	$(CC) -c -fopenmp BUSseq_alpha2.cpp -w -std=c++11

BUSseq_alpha3.o: BUSseq_alpha3.cpp
	$(CC) -c -fopenmp BUSseq_alpha3.cpp -w -std=c++11

BUSseq_alpha4.o: BUSseq_alpha4.cpp
	$(CC) -c -fopenmp BUSseq_alpha4.cpp -w -std=c++11

# clean all files
clean:
	rm -f *.o;