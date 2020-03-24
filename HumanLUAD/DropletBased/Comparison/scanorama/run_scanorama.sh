#!/bin/sh
cd DropletBased/Comparison/scanorama

# Conduct analysis by Scanorama
# Activate conda env
source activate tensorflow_gpuenv
# Python version 3.7.3
python scanorama_log_LUAD_v1.py
R --vanilla --slave < scanorama_log_LUAD_v1.R