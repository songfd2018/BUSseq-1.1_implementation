#!/bin/sh
cd DropletBased

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_KX" and the posterior inference in the folder "Inference_KX"
../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p LUAD -v 1 -K 2 -i 8000 -o 2000 -s 4629 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p LUAD -v 1 -K 2 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p LUAD -v 1 -K 3 -i 8000 -o 2000 -s 2548 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p LUAD -v 1 -K 3 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p LUAD -v 1 -K 4 -i 8000 -o 2000 -s 8880 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p LUAD -v 1 -K 4 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p LUAD -v 1 -K 5 -i 8000 -o 2000 -s 924 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p LUAD -v 1 -K 5 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p LUAD -v 1 -K 6 -i 8000 -o 2000 -s 798 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p LUAD -v 1 -K 6 -i 8000 -b 4000 -c 8

# Summarize the results of BUSseq
R --vanilla --slave < summarize_BUSseq.R
