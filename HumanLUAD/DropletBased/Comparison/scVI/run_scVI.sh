#!/bin/sh
cd DropletBased/Comparison/scVI

# Conduct analysis by scVI
# Activate conda environment
source activate tensorflow_env

# Python version 3.7.3
python scVI_LUAD_v1_log.py

