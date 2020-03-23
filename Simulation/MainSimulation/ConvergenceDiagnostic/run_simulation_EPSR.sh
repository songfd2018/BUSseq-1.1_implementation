#!/bin/sh
cd MainSimulation/ConvergenceDiagnostic

# run another BUSseq MCMC chain
../../../../BUSseq-1.0/BUSseq -d./ -r../RawCountData/ -p simulation -v 1 -K 5 -i 4000 -o 1000 -s 2187 -c 8
../../../../BUSseq-1.0/BUSseq_inference -d./ -r../RawCountData/ -p simulation -v 1 -K 5 -i 4000 -b 2000 -c 8

# Calculate EPSR
R --vanilla --slave < EPSR_analysis.R
