# normlization comparison
rm(list=ls())
library(DESeq2)
library(edgeR)
library(scran)
setwd("MainSimulation/Comparison")
# load the dimention of count data
dim <- unlist(read.table("../RawCountData/dim_simulation_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[4:(4+B-1)]

# load raw count data
rawcount <- as.matrix(read.table("../RawCountData/count_data_simulation_v1.txt"))
counts <- rawcount[,1:nb[1]]

# load the cell type indicator
w.est <- unlist(read.table("../True_para/w_syn.txt"))
cell_type_col <- c("black", "dodgerblue", "orange","chartreuse","blueviolet")
col <- cell_type_col[w.est[1:nb[1]] + 1]

dir.create("Images", showWarning=FALSE)

make.plot <- function(sf, truth, name, main="") {
  resids <- log2(sf) - log2(truth)
  fitted <- lm(resids ~ 1, na.action=na.exclude)
  all.range <- range(c(truth * 1.2,truth * 0.8))
  shuffle <- sample(1:length(sf))
  
  jpeg(paste0("Images/",name, ".jpg"), width = 600, height = 600)
  par(mar=c(5.1,5.1,4.1,1.1))
  plot(truth[shuffle], (truth * 2^residuals(fitted))[shuffle], xlim=all.range, ylim=all.range, 
       ylab="Estimated factors", xlab="True factors", log="xy", pch=1, cex = 2.5,
       col=col[shuffle], cex.axis=2.5, cex.lab=3, main=main, cex.main=3)
  abline(0, 1, col="red",lty = 2)
  dev.off()
}

# true values
factor_true <- exp(unlist(read.table("../True_para/delta_syn.txt"))) 
factor_true <- factor_true[1:nb[1]]

# BUSseq
delta.est <- unlist(read.table("../Inference_K5/delta_est.txt"))
delta.est <- delta.est[1:nb[1]]
make.plot(exp(delta.est), factor_true, "BUSseq_factor_est", main="BUSseq")

# library size
lib.sf <- colSums(counts)
make.plot(lib.sf, factor_true, "library_factor_est", main="Library size")

# DESeq
logvals <- log(counts)
logvals[is.infinite(logvals)] <- NA_real_
gm <- exp(rowMeans(logvals, na.rm=TRUE))
size.sf <- estimateSizeFactorsForMatrix(counts, geoMeans=gm)
make.plot(size.sf, factor_true, "DESseq_factor_est", main="DESeq")

# TMM
tmm.sf <- calcNormFactors(counts) * colSums(counts)
make.plot(tmm.sf, factor_true, "tmm_factor_est", main="TMM")

# deconvolution method
# Size factors with summation:
final.sf <- computeSumFactors(counts, clusters=NULL, min.mean=0)
make.plot(final.sf, factor_true, "Deconvolution_factor_est", main="Deconvolution")

# Size factors with clustering prior to summation:
emp.clusters <- quickCluster(counts)
final2.sf <- computeSumFactors(counts, clusters=emp.clusters, min.mean=0)
make.plot(final2.sf, factor_true, "Deconvolution_cluster_factor_est", main="Deconvolution with clustering")
