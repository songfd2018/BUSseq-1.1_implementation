#!/bin/sh

cd OtherComparison/FilterOutTEG

# filter out the common TEG
R --vanilla --slave < Hemat_TEG_filtering.R

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_K6" and the posterior inference in the folder "Inference_K6" 

../../../../BUSseq-1.0/BUSseq -d./  -r./RawCountData/ -p hemat -v 3 -K 6 -i 8000 -o 2000 -s 372 -c 8
../../../../BUSseq-1.0/BUSseq_inference -d./  -r./RawCountData/ -p hemat -v 3 -K 6 -i 8000 -b 4000 -c 8

# compare the clustering results
R --vanilla --slave < Gene_filter_comparison_hemat.R