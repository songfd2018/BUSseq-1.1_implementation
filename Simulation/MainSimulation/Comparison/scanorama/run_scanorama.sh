#!/bin/sh

# Conduct analysis by scanorama

cd MainSimulation/Comparison/scanorama
# Activate conda env
source activate tensorflow_gpuenv
# Python version 3.7.3
python scanorama_log_v1.py
R --vanilla --slave < scanorama_log_v1.R

# Summarize the results
R --vanilla --slave < summarize_scanorama.R