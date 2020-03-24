#!/bin/sh
cd DifferentSampleComplexity


# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_KX" and the posterior inference in the folder "Inference_KX"
../../../BUSseq-1.0/BUSseq -d./three_mix/ -r./RawCountData/ -p LUAD -v 2 -K 2 -i 8000 -o 2000 -s 7700 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./three_mix/ -r./RawCountData/ -p LUAD -v 2 -K 2 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./three_mix/ -r./RawCountData/ -p LUAD -v 2 -K 3 -i 8000 -o 2000 -s 7800 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./three_mix/ -r./RawCountData/ -p LUAD -v 2 -K 3 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./three_mix/ -r./RawCountData/ -p LUAD -v 2 -K 4 -i 8000 -o 2000 -s 2657 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./three_mix/ -r./RawCountData/ -p LUAD -v 2 -K 4 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./three_mix/ -r./RawCountData/ -p LUAD -v 2 -K 5 -i 8000 -o 2000 -s 3400 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./three_mix/ -r./RawCountData/ -p LUAD -v 2 -K 5 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./three_mix/ -r./RawCountData/ -p LUAD -v 2 -K 6 -i 8000 -o 2000 -s 2229 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./three_mix/ -r./RawCountData/ -p LUAD -v 2 -K 6 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./five_mix/ -r./RawCountData/ -p LUAD -v 3 -K 2 -i 8000 -o 2000 -s 7700 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./five_mix/ -r./RawCountData/ -p LUAD -v 3 -K 2 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./five_mix/ -r./RawCountData/ -p LUAD -v 3 -K 3 -i 8000 -o 2000 -s 1359 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./five_mix/ -r./RawCountData/ -p LUAD -v 3 -K 3 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./five_mix/ -r./RawCountData/ -p LUAD -v 3 -K 4 -i 8000 -o 2000 -s 7744 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./five_mix/ -r./RawCountData/ -p LUAD -v 3 -K 4 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./five_mix/ -r./RawCountData/ -p LUAD -v 3 -K 5 -i 8000 -o 2000 -s 3400 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./five_mix/ -r./RawCountData/ -p LUAD -v 3 -K 5 -i 8000 -b 4000 -c 8

../../../BUSseq-1.0/BUSseq -d./five_mix/ -r./RawCountData/ -p LUAD -v 3 -K 6 -i 8000 -o 2000 -s 8823 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./five_mix/ -r./RawCountData/ -p LUAD -v 3 -K 6 -i 8000 -b 4000 -c 8
