CC=g++

all: BUSseq_delta1 BUSseq_delta2 BUSseq_delta3 BUSseq_delta4

# link
BUSseq_delta1: BUSseq_delta1.o rngstream.o
	$(CC) -fopenmp BUSseq_delta1.o rngstream.o -o BUSseq_delta1

BUSseq_delta2: BUSseq_delta2.o rngstream.o
	$(CC) -fopenmp BUSseq_delta2.o rngstream.o -o BUSseq_delta2

BUSseq_delta3: BUSseq_delta3.o rngstream.o
	$(CC) -fopenmp BUSseq_delta3.o rngstream.o -o BUSseq_delta3

BUSseq_delta4: BUSseq_delta4.o rngstream.o
	$(CC) -fopenmp BUSseq_delta4.o rngstream.o -o BUSseq_delta4



# compile
BUSseq_delta1.o: BUSseq_delta1.cpp
	$(CC) -c -fopenmp BUSseq_delta1.cpp -w -std=c++11

BUSseq_delta2.o: BUSseq_delta2.cpp
	$(CC) -c -fopenmp BUSseq_delta2.cpp -w -std=c++11

BUSseq_delta3.o: BUSseq_delta3.cpp
	$(CC) -c -fopenmp BUSseq_delta3.cpp -w -std=c++11

BUSseq_delta4.o: BUSseq_delta4.cpp
	$(CC) -c -fopenmp BUSseq_delta4.cpp -w -std=c++11

# clean all files
clean:
	rm -f *.o;
