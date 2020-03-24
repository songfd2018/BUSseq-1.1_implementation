#!/bin/sh
cd ConvergenceDiagnostic

# run another BUSseq MCMC chain
../../../BUSseq-1.0/BUSseq -d./ -r../RawCountData/ -p hemat -v 1 -K 6 -i 8000 -o 2000 -s 372 -c 8
../../../BUSseq-1.0/BUSseq_inference -d./ -r../RawCountData/ -p hemat -v 1 -K 6 -i 8000 -b 4000 -c 8

# Calculate EPSR
R --vanilla --slave < EPSR_analysis.R
