#!/bin/sh
cd MainSimulation

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_KX", where X = 3 - 10, and the posterior inference in the folder "Inference_KX" 

../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 1 -K 3 -i 4000 -o 1000 -s 9722 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 1 -K 3 -i 4000 -b 2000 -c 8

../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 1 -K 4 -i 4000 -o 1000 -s 2088 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 1 -K 4 -i 4000 -b 2000 -c 8

../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 1 -K 5 -i 4000 -o 1000 -s 9123 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 1 -K 5 -i 4000 -b 2000 -c 8

../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 1 -K 6 -i 4000 -o 1000 -s 9641 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 1 -K 6 -i 4000 -b 2000 -c 8

../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 1 -K 7 -i 4000 -o 1000 -s 4618 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 1 -K 7 -i 4000 -b 2000 -c 8

../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 1 -K 8 -i 4000 -o 1000 -s 7245 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 1 -K 8 -i 4000 -b 2000 -c 8

../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 1 -K 9 -i 4000 -o 1000 -s 8200 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 1 -K 9 -i 4000 -b 2000 -c 8

../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 1 -K 10 -i 4000 -o 1000 -s 4191 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 1 -K 10 -i 4000 -b 2000 -c 8

