#!/bin/sh
cd RareCellType

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_KX" and the posterior inference in the folder "Inference_KX"
../../../BUSseq-1.0/BUSseq -d./v4/ -r./RawCountData/ -p LUAD -v 4 -K 3 -i 8000 -o 2000 -s 6252 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./v4/ -r./RawCountData/ -p LUAD -v 4 -K 3 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./v5/ -r./RawCountData/ -p LUAD -v 5 -K 3 -i 8000 -o 2000 -s 6252 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./v5/ -r./RawCountData/ -p LUAD -v 5 -K 3 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./v6/ -r./RawCountData/ -p LUAD -v 6 -K 3 -i 8000 -o 2000 -s 9394 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./v6/ -r./RawCountData/ -p LUAD -v 6 -K 3 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./v7/ -r./RawCountData/ -p LUAD -v 7 -K 3 -i 8000 -o 2000 -s 9394 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./v7/ -r./RawCountData/ -p LUAD -v 7 -K 3 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./v8/ -r./RawCountData/ -p LUAD -v 8 -K 3 -i 8000 -o 2000 -s 3064 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./v8/ -r./RawCountData/ -p LUAD -v 8 -K 3 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./v9/ -r./RawCountData/ -p LUAD -v 9 -K 3 -i 8000 -o 2000 -s 9394 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./v9/ -r./RawCountData/ -p LUAD -v 9 -K 3 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./v10/ -r./RawCountData/ -p LUAD -v 10 -K 3 -i 8000 -o 2000 -s 3064 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./v10/ -r./RawCountData/ -p LUAD -v 10 -K 3 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./v11/ -r./RawCountData/ -p LUAD -v 11 -K 3 -i 8000 -o 2000 -s 9394 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./v11/ -r./RawCountData/ -p LUAD -v 11 -K 3 -i 8000 -b 4000 -c 8
