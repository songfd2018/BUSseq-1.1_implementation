# Generate Figures for the Pancreas Study

This folder contains scripts that generate figures for the pancreas study in the manuscript. 

## BUSseq

First, please download the raw count data from [the OneDrive link](https://mycuhk-my.sharepoint.com/:u:/g/personal/1155082896_link_cuhk_edu_hk/EW-UIYqDLIRNk8DwWP823AUBKmeh_c9__Rs-7YrtOO34zA?e=PJlW1u) with password `BUSseq2019`. Then, please run `run_BUSseq.sh` to run BUSseq, output the posterior samples in the folder `MCMC_sampling_K8`, and conduct statistical inference in the folder `Inference_K8`. Next, please run `BUSseq/summarize_BUSseq.R` to correct count data and draw t-SNE and PCA plots in our manuscript.

We evaluate the fitness of BUSseq model to the observed data by calculating zero rates and dropout rates of all batches, conducting posterior predictive check on batch-specific zero rates and drawing mean-variance trends of all genes.

### Dropout rate

BUSseq can detect the existence of dropout events automatically. To calculate dropout rates in the pancreas study, please run `BUSseq/Dropout_rate_pancreas.R` to output the batch-specific observed zero rates and estimated dropout rates.

### Posterior predictive check 

To compare BUSseq with its simplified form (BUSseq with no zero inflation, `BUSseq_nzf`), we conduct posterior prediction check on zero rates to check which model fits the observed scRNA-seq data better. More specifically, we ran 8,000 iterations with the first 4,000 iterations as burn-ins in the MCMC algorithm, so we generated *J=8,000-4,000=4,000* replicated datasets. For each generated replicate dataset, we calculated the zero rate of each batch. Finally, we averaged the zero rates over all *J* iterations to calculate the posterior mean of the zero rate of each batch and compared it with the corresponding observed zero rate.

Please run `BUSseq/PPC_pancreas_original_zerorate.R` to generate the zero rate of each batch for all *J* iterations in the text file `BUSseq/ZeroRate_pancreas_v1.txt`. To make a comparison with `BUSseq_nzf`, please refer to the [BUSseq_nzf](#busseq_nzf) section.


### Mean-variance trend

In order to explore how well a model recapitulates the properties of real scRNA-seq data, Zappia et al. [[1]](#1) suggested simulating a scRNA-seq dataset and comparing its properties with those of real data. As batch one is taken as the reference batch in our model, without loss of generality, here we draw the mean-variance trends for the second batch *b=2* of the pancreas study as illustrations. Please run `BUSseq/Mean_variance_trend_comparison.R` to generate the mean-variance trends across all cell types and within each cell type.

## Benchmarked batch-effects-correction method

Each benchmarked method has its own clustering approach, as outlined in its respective publication. We first respected the original clustering algorithm and analysis protocol of each benchmarked method in our main manuscript. For comparison, we have also applied robust k-means clustering by the `pam` function of the R package [cluster](https://github.com/cran/cluster) as the clustering method for all the benchmarked batch-effects-correction methods. To calculate ARI and draw t-SNE and PCA plots of corrected data by each method, please

   - Run `BatchEffectCorrectionComparison/liger/run_liger.R` to apply [LIGER](https://github.com/MacoskoLab/liger) to the pancreas dataset.

   - Run `BatchEffectCorrectionComparison/MNN/run_MNN.R` to apply [MNN](https://github.com/MarioniLab/MNN2017) to the pancreas dataset.

   - Run `BatchEffectCorrectionComparison/scanorama/run_scanorama.sh` to apply [Scanorama](https://github.com/brianhie/scanorama) to the pancreas dataset.

   - Run `BatchEffectCorrectionComparison/scVI/run_scVI.sh` to apply [scVI](https://github.com/YosefLab/scVI) to the pancreas dataset.

   - Run `BatchEffectCorrectionComparison/Seurat/run_Seurat.R` to apply [Seurat](https://satijalab.org/seurat/) to the pancreas dataset.

   - Run `BatchEffectCorrectionComparison/ZINBWaVE/run_ZINBWaVE.R` to apply [ZINBWaVE](https://github.com/drisso/zinbwave) to the pancreas dataset.

Finally, run `collect_evaluation.R` to collect the ARIs of all of the methods and calculate Silhouette coefficients based on the t-SNE coordinates of each compared method. Meanwhile, t-SNE and PCA plots of uncorrected count data and the boxplot of Silhouette coefficients are also generated.

## Convergence diagnostic

We diagnose the convergence of BUSseq by acceptance rates of Metropolis-Hastings (MH) updates and estimated potential scale reduction (EPSR) factors.

### Acceptance rate

To calculate the acceptance rates of the MCMC algorithm, please run `ConvergenceDiagnostic/Acceptance_rate_cal_pancreas.R`. The code will output the acceptance rates of all parameters updated by MH steps.

### EPSR

To calculate EPSR of BUSseq, please run `ConvergenceDiagnostic/run_pancreas_EPSR.sh`. It first starts another chain and then calculates EPSR factors of these two chains by `ConvergenceDiagnostic/EPSR_analysis.R`. The R code will output the proportion of EPSRs less than 1.3 for mean expression levels of each cell type, location batch effects and overdispersion parameters. If all these three proportions are larger than 80\%, then the MCMC algorithm is regarded to have reached convergence.


## Others

### BUSseq_nzf

Here, we consider the simplified BUSseq model with no zero inflation, abbreviated as `BUSseq_nzf`. We compare `BUSseq` and `BUSseq_nzf` in terms of both model fitting, using a posterior predictive check for zero rates, and clustering accuracy with ARIs.

   - First, please run `OtherComparison/BUSseq_nzf/run_BUSseq_nzf.sh` to compile the C++ source code for `BUSseq_nzf` model and refit the pancreas dataset by `BUSseq_nzf` model. 

   - Then, please run `OtherComparison/BUSseq_nzf/PPC_pancreas_nodropout_zerorate.R` to generate the zero rates of batches for all *J* iterations in the txt file `OtherComparison/BUSseq_nzf/ZeroRate_pancreas_nzf.txt`. 
   - Finally, please run `OtherComparison/BUSseq_nzf/Model_Comparison.R` to compare the posterior predictive checks on zero rates and ARIs between `BUSseq` and `BUSseq_nzf`.

### Filtering genes according to expression levels

In addition to filtering highly variable genes (HVG), we also followed Duò et al. [[2]](#2) in selecting the genes with the highest mean expression levels across all cells. We retained genes with the top 50\% highest average of log-scale expression values such that 408 common genes were retained from the four pancreas batches for downstream analysis. We denote these genes as top ranked genes in terms of expression levels (TEG).

Please first download the raw count data of pancreas study after quality control from [the OneDrive link](https://mycuhk-my.sharepoint.com/:u:/g/personal/1155082896_link_cuhk_edu_hk/EbNh9zPYMDhKiOKNLJo9azwBBlhUpjMwbjwhfLCwPKI5Ww?e=6nhveg) with password `BUSseq2020` and unzip it to `OtherComparison/FilterOutTEG/Data`. Then, please run `OtherComparison/FilterOutTEG/run_BUSseq_HEG.sh` directly. 

   - First, `OtherComparison/FilterOutTEG/Pancreas_TEG_filtering.R` filters out common TEGs across four batches.

   - Then, BUSseq is applied to these TEGs.

   - Finally, `OtherComparison/FilterOutTEG/Gene_filter_comparison_pancreas.R` compares cell type labeling between those output by BUSseq when applied to HVGs and those when applied to TEGs.

## References
<a id="1">[1]</a> 
Zappia, L., Phipson, B. and Oshlack, A., 2017. Splatter: simulation of single-cell RNA sequencing data. Genome Biology, 18(1), p.174.
<a id="2">[2]</a>
Duò, A., Robinson, M.D. and Soneson, C., 2018. A systematic performance evaluation of clustering methods for single-cell RNA-seq data. F1000Research, 7.
