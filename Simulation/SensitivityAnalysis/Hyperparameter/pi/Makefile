CC=g++

all: BUSseq_pi1 BUSseq_pi2 BUSseq_pi3 BUSseq_pi4

# link
BUSseq_pi1: BUSseq_pi1.o rngstream.o
	$(CC) -fopenmp BUSseq_pi1.o rngstream.o -o BUSseq_pi1

BUSseq_pi2: BUSseq_pi2.o rngstream.o
	$(CC) -fopenmp BUSseq_pi2.o rngstream.o -o BUSseq_pi2

BUSseq_pi3: BUSseq_pi3.o rngstream.o
	$(CC) -fopenmp BUSseq_pi3.o rngstream.o -o BUSseq_pi3

BUSseq_pi4: BUSseq_pi4.o rngstream.o
	$(CC) -fopenmp BUSseq_pi4.o rngstream.o -o BUSseq_pi4



# compile
BUSseq_pi1.o: BUSseq_pi1.cpp
	$(CC) -c -fopenmp BUSseq_pi1.cpp -w -std=c++11

BUSseq_pi2.o: BUSseq_pi2.cpp
	$(CC) -c -fopenmp BUSseq_pi2.cpp -w -std=c++11

BUSseq_pi3.o: BUSseq_pi3.cpp
	$(CC) -c -fopenmp BUSseq_pi3.cpp -w -std=c++11

BUSseq_pi4.o: BUSseq_pi4.cpp
	$(CC) -c -fopenmp BUSseq_pi4.cpp -w -std=c++11

# clean all files
clean:
	rm -f *.o;