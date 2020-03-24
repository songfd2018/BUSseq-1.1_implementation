#!/bin/sh
cd BatchEffectCorrectionComparison/scanorama
# Conduct analysis by Scanorama
# Activate conda env
source activate tensorflow_gpuenv
# Python version 3.7.3
python scanorama_log.py
R --vanilla --slave < scanorama_log.R

# Summarize the results
R --vanilla --slave < summarize_scanorama.R