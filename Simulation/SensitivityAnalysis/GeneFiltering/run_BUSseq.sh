#!/bin/sh
cd SensitivityAnalysis/GeneFiltering

# simulate data
R --vanilla --slave < simulate_data_gene_filtering.R

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_KX" and the posterior inference in the folder "Inference_KX" 
../../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 3 -K 3 -i 4000 -o 1000 -s 5101 -c 8
../../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 3 -K 3 -i 4000 -b 2000 -c 8

../../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 3 -K 4 -i 4000 -o 1000 -s 4473 -c 8
../../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 3 -K 4 -i 4000 -b 2000 -c 8

../../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 3 -K 5 -i 4000 -o 1000 -s 9123 -c 8
../../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 3 -K 5 -i 4000 -b 2000 -c 8

../../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 3 -K 6 -i 4000 -o 1000 -s 6390 -c 8
../../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 3 -K 6 -i 4000 -b 2000 -c 8

../../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 3 -K 7 -i 4000 -o 1000 -s 5614 -c 8
../../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 3 -K 7 -i 4000 -b 2000 -c 8

../../../../BUSseq-1.0/BUSseq -d./ -r./RawCountData/ -p simulation -v 3 -K 8 -i 4000 -o 1000 -s 7245 -c 8
../../../../BUSseq-1.0/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 3 -K 8 -i 4000 -b 2000 -c 8

# Summarize the results of BUSseq
R --vanilla --slave < summarize_BUSseq.R
