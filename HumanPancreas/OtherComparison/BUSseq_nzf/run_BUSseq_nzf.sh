#!/bin/sh
cd OtherComparison/BUSseq_nzf

# compile the source code
make

# run BUSseq with no zero inflation 
BUSseq_nd -d./  -r../../RawCountData/ -p pancreas -v 1 -K 8 -i 8000 -o 2000 -s 4537 -c 8
BUSseq_inference_nd -d./  -r../../RawCountData/ -p pancreas -v 1 -K 8 -i 8000 -b 4000 -c 8
