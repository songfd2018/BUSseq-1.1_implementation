# Analyze the results of the posterior predictive check
rm(list=ls())
library(mclust) # For ARI

setwd("OtherComparison/BUSseq_nzf")
proj <- "pancreas"

# load the raw count data
y_obs <- read.table(paste0("../../RawCountData/count_data_",proj,"_v1.txt"))

# load dimension information
dim <- unlist(read.table(paste0("../../RawCountData/dim_",proj,"_v1.txt")))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table(paste0("../../RawCountData/metadata_",proj,"_v1.txt"))

# calculate the observed zero rate
zerorate_obs <- sum(y_obs == 0) / N / G

# calculate the observed zero rate within each batch
zerorate_batch <- rep(NA,B)
ind_first <- 0
for(b in 1:B){
  ind_batch <- ind_first + 1:nb[b]
  zerorate_batch[b] <- sum(y_obs[,ind_batch] == 0) / nb[b] / G
  ind_first <- ind_first + nb[b]
}

# load the zero rate of no dropout
zerorate_rep_nd <- read.table("ZeroRate_pancreas_nzf.txt")
summary_nd <- matrix(NA, 6, B + 1)
for(i in 1:(B+1)){
  summary_nd[,i] <- c(quantile(zerorate_rep_nd[,i],probs = c(0,0.25,0.5)),
                      mean(zerorate_rep_nd[,i]), 
                      quantile(zerorate_rep_nd[,i], probs = c(0.75, 1)))
}
print(rbind(c(zerorate_obs, zerorate_batch),summary_nd))


# load the zero rate with dropout
zerorate_rep_or <- read.table("../../BUSseq/ZeroRate_pancreas_v1.txt")
summary_or <- matrix(NA, 6, B + 1)
for(i in 1:(B+1)){
  summary_or[,i] <- c(quantile(zerorate_rep_or[,i],probs = c(0,0.25,0.5)),
                      mean(zerorate_rep_or[,i]), 
                      quantile(zerorate_rep_or[,i], probs = c(0.75, 1)))
}
print(rbind(c(zerorate_obs, zerorate_batch),summary_or))

# ARI comparison
w_BUSseq <- unlist(read.table("../../Inference_K8/w_est.txt"))
ARI_BUSseq <- adjustedRandIndex(metadata$CellType,w_BUSseq)

w_BUSseq_nzf <- unlist(read.table("Inference_K8/w_est.txt"))
ARI_BUSseq_nzf <- adjustedRandIndex(metadata$CellType,w_BUSseq_nzf)

message(paste0("The ARI of BUSseq is ",ARI_BUSseq,", while the ARI of BUSseq_nzf is ",ARI_BUSseq_nzf,"."))

