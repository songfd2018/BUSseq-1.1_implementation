#!/bin/sh

cd OtherComparison/FilterOutTEG

# conduct gene filtering
R --vanilla --slave < Pancreas_HEG_filtering.R

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_K5" and the posterior inference in the folder "Inference_K5" 

../../../../BUSseq-1.0/BUSseq -d./  -r./RawCountData/ -p pancreas -v 3 -K 8 -i 8000 -o 2000 -s 3268 -c 8
../../../../BUSseq-1.0/BUSseq_inference -d./  -r./RawCountData/ -p pancreas -v 3 -K 8 -i 8000 -b 4000 -c 8

# compare the clustering results
R --vanilla --slave < Gene_filter_comparison_pancreas.R