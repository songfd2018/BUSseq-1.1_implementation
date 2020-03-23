CC=g++

all: BUSseq_p1 BUSseq_p2 BUSseq_p3 BUSseq_p4

# link
BUSseq_p1: BUSseq_p1.o rngstream.o
	$(CC) -fopenmp BUSseq_p1.o rngstream.o -o BUSseq_p1

BUSseq_p2: BUSseq_p2.o rngstream.o
	$(CC) -fopenmp BUSseq_p2.o rngstream.o -o BUSseq_p2

BUSseq_p3: BUSseq_p3.o rngstream.o
	$(CC) -fopenmp BUSseq_p3.o rngstream.o -o BUSseq_p3

BUSseq_p4: BUSseq_p4.o rngstream.o
	$(CC) -fopenmp BUSseq_p4.o rngstream.o -o BUSseq_p4



# compile
BUSseq_p1.o: BUSseq_p1.cpp
	$(CC) -c -fopenmp BUSseq_p1.cpp -w -std=c++11

BUSseq_p2.o: BUSseq_p2.cpp
	$(CC) -c -fopenmp BUSseq_p2.cpp -w -std=c++11

BUSseq_p3.o: BUSseq_p3.cpp
	$(CC) -c -fopenmp BUSseq_p3.cpp -w -std=c++11

BUSseq_p4.o: BUSseq_p4.cpp
	$(CC) -c -fopenmp BUSseq_p4.cpp -w -std=c++11

# clean all files
clean:
	rm -f *.o;
