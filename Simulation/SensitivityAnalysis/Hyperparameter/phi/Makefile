CC=g++

all: BUSseq_phi1 BUSseq_phi2 BUSseq_phi3 BUSseq_phi4

# link
BUSseq_phi1: BUSseq_phi1.o rngstream.o
	$(CC) -fopenmp BUSseq_phi1.o rngstream.o -o BUSseq_phi1

BUSseq_phi2: BUSseq_phi2.o rngstream.o
	$(CC) -fopenmp BUSseq_phi2.o rngstream.o -o BUSseq_phi2

BUSseq_phi3: BUSseq_phi3.o rngstream.o
	$(CC) -fopenmp BUSseq_phi3.o rngstream.o -o BUSseq_phi3

BUSseq_phi4: BUSseq_phi4.o rngstream.o
	$(CC) -fopenmp BUSseq_phi4.o rngstream.o -o BUSseq_phi4



# compile
BUSseq_phi1.o: BUSseq_phi1.cpp
	$(CC) -c -fopenmp BUSseq_phi1.cpp -w -std=c++11

BUSseq_phi2.o: BUSseq_phi2.cpp
	$(CC) -c -fopenmp BUSseq_phi2.cpp -w -std=c++11

BUSseq_phi3.o: BUSseq_phi3.cpp
	$(CC) -c -fopenmp BUSseq_phi3.cpp -w -std=c++11

BUSseq_phi4.o: BUSseq_phi4.cpp
	$(CC) -c -fopenmp BUSseq_phi4.cpp -w -std=c++11

# clean all files
clean:
	rm -f *.o;
