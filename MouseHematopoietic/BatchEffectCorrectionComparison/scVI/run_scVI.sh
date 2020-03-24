#!/bin/sh
cd BatchEffectCorrectionComparison/scVI
# Conduct analysis by scVI
# Activate conda environment
source activate tensorflow_env
# Python version 3.7.3
python scVI_log.py
# Summarize the results
R-351 --vanilla < summarize_scVI.R
