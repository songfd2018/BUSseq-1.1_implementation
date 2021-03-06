﻿# Generate Figures for the Simulation Studies

This folder contains scripts that generate figures for the simulation studies in the manuscript.

## Main simulation

To generate the main simulation data and apply BUSseq on it, please:

   - Run `MainSimulation/simulate_data.R` to generate a synthetic dataset together with dimension information and metadata in the `RawCountData` folder

   - Run `MainSimulation/run_BUSseq.sh` to obtain the posterior samples and conduct statistical inference.

### Convergence diagnostic

We diagnose the convergence of BUSseq by acceptance rates of Metropolis-Hastings (MH) updates and estimated potential scale reduction (EPSR) factors.

To calculate the acceptance rate of MCMC algorithm, please run `MainSimulation/ConvergenceDiagnostic/Acceptance_rate_cal_simulation.R`. The code will output the acceptance rate of all parameters updated by MH steps.

To calculate EPSR of BUSseq, please run `MainSimulation/ConvergenceDiagnostic/run_simulation_EPSR.sh`. It first starts another chain and then calculates EPSR factors of these two chains by `MainSimulation/ConvergenceDiagnostic/EPSR_analysis.R`. The R codes will output the proportion of EPSRs less than 1.3 for mean expression levels, location batch effects and overdispersion parameters. If all these three proportions are larger than 80\%, then the MCMC algorithm is regarded to have reached convergence.

### Comparison

First, please run `MainSimulation/BUSseq/summarize_BUSseq.R` to obtain the parameter estimates and the corrected count data.

#### Batch effect correction

Each benchmarked method has its own clustering approach, as outlined in its respective publication. We respect the original clustering algorithm and analysis protocol of each benchmarked method. To calculate ARI and draw t-SNE and PCA plots of corrected data by each method, please

   - Run `MainSimulation/Comparison/liger/run_liger.R` to apply [LIGER](https://github.com/MacoskoLab/liger) to the simulation dataset.

   - Run `MainSimulation/Comparison/MNN/run_MNN.R` to apply [MNN](https://github.com/MarioniLab/MNN2017) to the simulation dataset.

   - Run `MainSimulation/Comparison/scanorama/run_scanorama.sh` to apply [Scanorama](https://github.com/brianhie/scanorama) to the simulation dataset.

   - Run `MainSimulation/Comparison/scVI/run_scVI.sh` to apply [scVI](https://github.com/YosefLab/scVI) to the simulation dataset.

   - Run `MainSimulation/Comparison/Seurat/run_Seurat.R` to apply [Seurat](https://satijalab.org/seurat/) to the simulation dataset.

   - Run `MainSimulation/Comparison/ZINBWaVE/run_ZINBWaVE.R` to apply [ZINBWaVE](https://github.com/drisso/zinbwave) to the simulation dataset.

Finally, please run `MainSimulation/collect_evaluation.R` to collect the ARIs of all methods and calculate Silhouette coefficients based on the t-SNE coordinates of each compared method. Meanwhile, t-SNE and PCA plots of uncorrected count data and the boxplot of Silhouette coefficients are generated.

#### Normalization

To compare BUSseq with the state-of-the-art normalization methods, including [DESeq](https://github.com/mikelove/DESeq2) normalization, [trimmed mean of M-values](https://github.com/StoreyLab/edge) normalization, [deconvolution](https://github.com/Albluca/scran) normalization, library size normalization, we draw scatter plots of the estimates of cell size factors output by each method against the true cell size factors. The true values of size factors are known in the simulation dataset, so we apply all these methods to the simulation dataset. Moreover, to avoid the impact of batch effects, we focused on the first batch with b = 1. Please run `MainSimulation/Comparison/normlization_comparison.R` to generate scatter plots in the folder `MainSimulation/Comparison/Images`.

#### Imputation

To compare BUSseq with the state-of-the-art imputation methods, including [SAVER](https://github.com/mohuangx/SAVER), [DrImpute](https://github.com/gongx030/DrImpute) and [scImpute](https://github.com/Vivianstats/scImpute), we calculate the Euclidean distance between the imputed values and the underlying true read counts for all of the observed zeros. Please run `MainSimulation/Comparison/impute_comparison.R` to output the Euclidean distance and zero rates of all methods.

## Sensitivity Analysis

### BUSseq is robust to the choice of hyper-parameters

For each hyperparameter, we vary its value to four different levels while fixing the other hyperparameters. Please run `SensitivityAnalysis/Hyperparameter/XXX/run_BUSseq_XXX.sh` to apply BUSseq with different values of a given hyperparameter, where `XXX` denotes the corresponding parameter whose prior distribution involves this hyperparameter. `Collect_inference_XXX.R` will output the four ARIs. If ARI is equal to one, then BUSseq perfectly cluster cells by cell type.

### BUSseq is robust to high zero rates

Here, we construct a simulation dataset with high zero rates. More specifically, the zero rates of the four batches are 53.55\%, 61.13\%, 60.95\% and 80.99\%, respectively. To generate the figures in the manuscript, please run `SensitivityAnalysis/HighZeroRate/run_BUSseq.sh` to:

   - generate the simulation dataset;

   - draw heatmaps of true parameters;

   - apply BUSseq to the simulation dataset;

   - draw heatmaps of parameter estimates as well as t-SNE and PCA plots of the corrected count data.

### BUSseq is robust to model misspecification

We further apply BUSseq to a dataset where the model assumption is violated. In BUSseq, we assume that each gene's overdispersion parameter is the same for all of the cells measured on a given batch; in the new simulation setting, we let the overdispersion parameters vary across cell types. Specifically, we regard a gene as highly expressed in a given cell type if its expression level in that cell type is higher than its mean expression level across all cell types. To generate the figure in the manuscript, please run `SensitivityAnalysis/ModelMisspecification/run_BUSseq.sh` to:

   - generate the simulation dataset;

   - draw heatmaps of true parameters;

   - apply BUSseq to the simulation dataset;

   - draw heatmaps of parameter estimates as well as t-SNE and PCA plots of the corrected count data.

### BUSseq is robust to gene filtering

In real data analyses, researchers usually conduct gene filtering and focus on those highly variable genes (HVG). Therefore, to mimic the real practice, we simulated a dataset with *G=20,000* genes and a total number of *N = 1,000* cells of *K = 5* cell types, measured by *B = 4* batches. Among all *20,000* genes, there are *1,326* intrinsic genes, which are differentially expressed between at least two different cell types. After gene filtering, we obtained *3,461* common HVGs across four batches. To generate the figures in the manuscript, please directly run `SensitivityAnalysis/GeneFiltering/run_BUSseq.sh` to:

   - generate the simulation dataset;

   - draw heatmaps of true parameters;

   - apply BUSseq to the simulation dataset;
   
   - draw heatmaps of parameter estimates as well as t-SNE and PCA plots of the corrected count data.

## Scalability

We explore how the RAM consumption and the computing time vary with the number of genes *G*, the number of cells *N* and the number of iterations per hard-disk writing. Please run `Scalability\scatter_time_RAM.R` to regenerate the chart plots in our manuscript. The generated images will be stored in the folder `Scalability\Image`.



