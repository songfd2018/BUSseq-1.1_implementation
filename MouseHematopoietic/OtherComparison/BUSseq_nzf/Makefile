CC=g++

all: BUSseq_nd BUSseq_inference_nd

# link
BUSseq_nd: BUSseq_nd.o rngstream.o
	$(CC) -fopenmp BUSseq_nd.o rngstream.o -o BUSseq_nd

BUSseq_inference_nd: BUSseq_inference_nd.o rngstream.o
	$(CC) -fopenmp BUSseq_inference_nd.o rngstream.o -o BUSseq_inference_nd

# compile
BUSseq_nd.o: BUSseq_nd.cpp
	$(CC) -c -fopenmp BUSseq_nd.cpp -w -std=c++11

BUSseq_inference_nd.o: BUSseq_inference_nd.cpp
	$(CC) -c -fopenmp BUSseq_inference_nd.cpp -w -std=c++11

rngstream.o: rngstream.cpp
	$(CC) -c -fopenmp -lgomp rngstream.cpp

# clean all files
clean:
	rm -f *.o;