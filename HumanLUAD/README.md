# Generate Figures for the Lung Adenocarcinoma (LUAD) Study

This folder contains scripts that generate figures for the LUAD study in the manuscript. 

To obtain the raw count data, please first clone [`sc_mixology` GitHub repository](https://github.com/LuyiTian/sc_mixology) to the current directory.

## Application of BUSseq to droplet-based scRNA-seq datasets

We further analyzed a dataset that contains samples assayed by droplet-based scRNA-seq protocols. Tian et al. [[1]](#1) created scRNA-seq datasets with known cell type labels by profiling cells from cancer cell lines. In one experiment, they assayed three LUAD cell lines---HCC827, H1975 and H2228 on three platforms with CELseq2, 10x Chromium and Drop-seq protocols, respectively. The last two protocols are droplet-based. 

Please run `DropletBased/data_preprocessing.R` to conduct data preprocessing. The processed data will be stored in the folder `DropletBased/RawCountData`. Meanwhile, PCA and t-SNE plots of raw count data will be drawn.

### BUSseq

To apply BUSseq to the LUAD dataset, please directly run `DropletBased/run_BUSseq.sh` to:

   - Apply BUSseq and vary the number of cell types *K* from 2 to 6. The posterior samples and inference results are stored in the folder `MCMC_sampling_KX` and `Inference_KX`, respectively, X = 2,3,4,5,6.
   - Select the number of cell types by BIC, correct batch effects, calculate ARI and draw t-SNE and PCA plots as well as heatmap of mean expression levels by `DropletBased/summarize_BUSseq.R`.

### Comparison with batch-effects-correction methods

Each benchmarked method has its own clustering approach, as outlined in its respective publication. We respect the original clustering algorithm and analysis protocol of each benchmarked method. To calculate ARI of each method, please

   - Run `DropletBased/Comparison/scanorama/run_scanorama.sh` to apply [Scanorama](https://github.com/brianhie/scanorama) to the LUAD dataset and output the corrected data.

   - Run `DropletBased/Comparison/scVI/run_scVI.sh` to apply [scVI](https://github.com/YosefLab/scVI) to the LUAD dataset and output the corrected data.

   - Run `DropletBased/method_comparison.R` to apply [LIGER](https://github.com/MacoskoLab/liger), [MNN](https://github.com/MarioniLab/MNN2017, [Seurat](https://satijalab.org/seurat/) and [ZINBWaVE](https://github.com/drisso/zinbwave) to the LUAD dataset and analyze the corrected data by scanorama and scVI. It will output the ARIs of all the benchmarked methods.

### Comparison with CellAssign

Compared with all the other batch effect correction methods, [CellAssign](https://github.com/Irrationone/cellassign) requires knowledge of the number of cell types and a set of marker genes for each cell type. We still compare the performance of cell type labeling between BUSseq and CellAssign. For a fair comparison, instead of selecting the number of cell types by BIC, when we apply BUSseq to the data we also set the number of cell type as three, corresponding to HCC827, H1975 and H2228 cell lines.

Please run `DropletBased/CellAssign/Comparison_with_cellassign.R` to

   - Construct marker genes from bulk RNA-seq data following [the vigenette in CellAssign R package]{https://github.com/Irrationone/cellassign/blob/master/vignettes/constructing-markers-from-purified-data.Rmd}
   - Apply CellAssign on the selected marker genes to infer the cell type labeling of all cells
   - Compare the performance of labeling between BUSseq and CellAssign in terms of ARI, F1-score and accuracy 

## Performance of BIC under different levels of sample complexity

To assess the performance of selecting the number of cell types by BIC across different levels of sample complexity, we consider four batches of pseudo LUAD cells assayed by Tian et al.[[1]](#1). More specifically, in each batch, single cells from the three cell lines were sorted into 384-well plates with 9 cells per well in different combinations. RNAs from all of the 9 cells in the same well were then pooled and sub-sampled such that an approximately single-cell quantity of RNA was extracted from each well. A triad *(a,b,c)* to denote a combination of *a* HCC827 cells, *b* H1975 cells and *c* H2228 cells for each well. Pseudo cells are much more similar to each other than cells from the three pure cell lines, so we can assess the performance of BUSseq in the experiment with a complex cell population.

First, we applied BUSseq to three pure cell line mixtures, including (9,0,0), (0,9,0) and (0,0,9), in the four batches. Next, we incorporated two extra cell line mixtures, (4,0,5) and (5,0,4) and applied BUSseq to five cell line mixtures. To diagnose the performance of BIC on these two datasets, please run

   - `generate_three_mix_LUAD.R` and `generate_five_mix_LUAD.R` to conduct data preprocessing and draw t-SNE and PCA plots of raw count data on three pure cell line mixtures and five cell line mixtures, respectively.

   - `run_BUSseq_LUADmix.sh` to conduct statistical inference by BUSseq on these two datasets

   - `summarize_BUSseq.R` to draw the BIC plots, correct batch effects, calculate ARIs and draw t-SNE and PCA plots of corrected count data in our manuscript.


## Identification of rare cell type

We continued to consider the three pure cell lines, HCC827 for (9,0,0), H1975 for (0,9,0) and H2228 for (0,0,9) in the last section. To mimic unevenly distributed rare cell types, we subsampled HCC827 cells from each batch. To explore whether the *rare* cell type HCC827 can be identified by BUSseq, please run

   - `downsample_LUADmix.R` to generate raw count data after downsampling HCC827 cells for 7 different cell numbers.
   
   - `run_BUSseq.sh` to conduct statistical inference by BUSseq on all of these datasets.

   - `Collect_ARI.R` to collect their ARIs.

## References
<a id="1">[1]</a> Tian, L., Dong, X., Freytag, S., Lê Cao, K.A., Su, S., JalalAbadi, A., Amann-Zalcenstein, D., Weber, T.S., Seidi, A., Jabbari, J.S. and Naik, S.H., 2019. Benchmarking single cell RNA-sequencing analysis pipelines using mixture control experiments. Nature methods, 16(6), pp.479-487.