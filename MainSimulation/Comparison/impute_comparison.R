# imputation comparison
rm(list=ls())
setwd("MainSimulation/Comparison")
### This is a Wrapper for imputation functions to be used in other scripts. ###
# Impute_default_all is used in analyzing real and splatter simulated data.

# Format required to add additional imputation methods:
#	input: a SingleCellExperiment object with "counts" and "logcounts"
#		a parameter value (with default), 
#		whether the input should be normalized (with default=TRUE)
#		number of cores to use (if applicable)
#		random seed to ensure reproducibility (with default).
#	output: the imputed matrix where rows=genes, cols=cells

require("SingleCellExperiment")
require("DrImpute") 
require("scImpute") 
require("SAVER") 
require("scater")

scImpute_wrapper<- function(sce, param=0.5, do.norm=TRUE, seed=42) {
  # normalization is never used
  tmp_dir <- paste(paste("Tmp",round(runif(1)*100000), sep="_"),"/",sep="")
  dir.create(tmp_dir)
  set.seed(seed)
  param_name <- "Dropout Threshold"
  saveRDS(as.matrix(assays(sce)[["counts"]]), file=paste(tmp_dir,"tmp.rds", sep=""))
  scImpute::scimpute(paste(tmp_dir,"tmp.rds", sep=""), infile="rds", outfile="rds",
                     type="count", drop_thre=param, out_dir=tmp_dir,
                     Kcluster=4, ncores=1)
  out <- readRDS(paste(tmp_dir,"scimpute_count.rds", sep=""))
  unlink(tmp_dir, recursive=TRUE)
  return(out);
}

DrImpute_wrapper <- function(sce, param=0, do.norm=TRUE, seed=42) {
  # uses pre-defined log-normalized matrix
  set.seed(seed)
  param_name <- "Zeros Remaining"
  out <- DrImpute::DrImpute(as.matrix(assays(sce)[["logcounts"]]),
                            ks=4,
                            zerop=param)
  return(out)
}

SAVER_wrapper <- function(sce, param=1, do.norm=TRUE, seed=42, n.cores=16){
  # optional CPM-like normalization
  n.cores <- as.numeric(n.cores)
  set.seed(seed)
  require(doParallel)
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  param_name <- "Percent of Genes"
  sf <- 1
  if (do.norm) {
    sf <- NULL
  }
  if (param == 1) {
    out <- saver(assays(sce)[["counts"]], do.fast=TRUE, size.factor=sf, ncores=n.cores)
  } else {
    out <- saver(assays(sce)[["counts"]], do.fast=TRUE, size.factor=sf,
                 npred=nrow(sce)*param, ncores=n.cores)
  }
  stopCluster(cl)
  return(out$estimate);
}


set.seed(28198)
n.cores <- 4

# load the dimention of count data
dim <- unlist(read.table("../RawCountData/dim_simulation_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[4:(4+B-1)]

# load raw count data
rawcount <- as.matrix(read.table("../RawCountData/count_data_simulation_v1.txt"))
counts <- rawcount[,1:nb[1]]

rownames(counts) <- paste0("Gene_",1:G)
colnames(counts) <- paste0("Cell_",1:nb[1])

# conduct library size normalization
sf <- colSums(counts)
norm_counts <- t(t(counts)/sf*median(sf))
lognorm_counts <- log(norm_counts+1)


load("../simulation_countdata_v1.RData")
mu_syn <- alpha.syn + beta.syn
mean_level_syn <- apply(mu_syn, 1, mean)

underlying_counts <- x[[1]][,1:nb[1]]


sim_sce <- SingleCellExperiment(assays=list(counts=counts, underlying_counts = underlying_counts, logcounts=lognorm_counts))

dir.create("Images", showWarning=FALSE)

BUSseq_imputed <- as.matrix(read.table("../Inference_K5/x_imputed.txt"))
BUSseq_imputed <- BUSseq_imputed[,1:nb[1]]
rownames(BUSseq_imputed) <- paste0("Gene_",1:G)
colnames(BUSseq_imputed) <- paste0("Cell_",1:nb[1])
assays(sim_sce)[["busseq"]] <- BUSseq_imputed[,1:nb[1]]

res <- scImpute_wrapper(sim_sce, do.norm=FALSE)
assays(sim_sce)[["sci"]] <- res

res <- DrImpute_wrapper(sim_sce, do.norm=FALSE)
assays(sim_sce)[["dri"]] <- res

res <- SAVER_wrapper(sim_sce, n.cores=n.cores, do.norm=FALSE)
assays(sim_sce)[["saver"]] <- res

# round all the imputed data

BUSseq_imputed <- BUSseq_imputed[,1:nb[1]]

scImpute_imputed <- assays(sim_sce)[["sci"]]
scImpute_imputed_round <- round(scImpute_imputed)

DrImpute_imputed <- exp(assays(sim_sce)[["dri"]])-1
DrImpute_imputed <- t (t(DrImpute_imputed) * sf / median(sf))
DrImpute_imputed_round <- round(DrImpute_imputed)


SAVER_imputed <- t(t(round(assays(sim_sce)[["saver"]])) * sf / median(sf))
SAVER_imputed_round <- round(SAVER_imputed)

# zero rate 
mean(counts == 0)
mean(underlying_counts==0)
mean(BUSseq_imputed==0)
mean(scImpute_imputed_round==0)
mean(DrImpute_imputed_round==0)
mean(SAVER_imputed_round==0)

# Eudineal Distance
dist_obs <- sqrt(sum(log1p(counts[counts==0]) - log1p(underlying_counts[counts==0]))^2)
print(dist_obs )
dist_BUSseq <- sqrt(sum(log1p(BUSseq_imputed[counts==0]) - log1p(underlying_counts[counts==0]))^2)
print(dist_BUSseq)
dist_scImpute <- sqrt(sum(log1p(scImpute_imputed[counts==0]) - log1p(underlying_counts[counts==0]))^2)
print(dist_scImpute)
dist_DrImpute <- sqrt(sum(log1p(DrImpute_imputed[counts==0]) - log1p(underlying_counts[counts==0]))^2)
print(dist_DrImpute)
dist_SAVER <- sqrt(sum(log1p(SAVER_imputed[counts==0]) - log1p(underlying_counts[counts==0]))^2)
print(dist_SAVER)

save.image("Imputation_comparison.RData")