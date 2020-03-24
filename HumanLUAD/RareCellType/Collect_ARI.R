# collect the ARI of downsampling different numbers of HCC827 cells
setwd("RareCellType")
library(mclust)
proj <- "LUAD"

# known cell type number
K <- 3

ver <- 4:11
n_set <- length(ver)

ARI_values <- rep(NA, n_set)

for(i in 1:n_set){
  
  # Load the cell type labels
  metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver[i],".txt"))
  w_true <- factor(paste0(metadata$Num_Cell1,"-",metadata$Num_Cell2,"-",metadata$Num_Cell3))
  
  # Load the estimated cell type labels
  w_est <- unlist(read.table(paste0("v",ver[i],"/Inference_K",K,"/w_est.txt")))
  
  # record ARI  
  ARI_values[i] <- adjustedRandIndex(w_est,w_true)
}

print(ARI_values)