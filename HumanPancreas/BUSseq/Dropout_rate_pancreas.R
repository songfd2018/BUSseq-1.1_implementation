# calculate the dropout rate in each dataset.
rm(list=ls())

set.seed(12345)
proj <- "pancreas"
ver <- 1

#############
# Load Data #
#############
# Working directory
setwd("BUSseq")

# Loading pancreatic count data
y_obs <- read.table(paste0("../RawCountData/count_data_",proj,"_v",ver,".txt"))

# Load dimension
dim <- read.table(paste0("../RawCountData/dim_",proj,"_v",ver,".txt"))
dim <- unlist(dim)
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3+1:B]

# Load metadata
metadata <- read.table(paste0("../RawCountData/metadata_",proj,"_v",ver,".txt"))

# Load gene_list
gene_list <- unlist(read.table(paste0("../RawCountData/gene_list_",proj,"_v",ver,".txt"),stringsAsFactors = F))

# Load the imputed count data
x_imputed <- read.table("../Inference_K8/x_imputed.txt")

##############################################################
# Calculate the zero rate and dropout rate within each batch #
##############################################################
zero_rate <- rep(NA, B)

drop_rate <- rep(NA, B)

ind_first <- 0
for(b in 1:B){
  ind_batch <- ind_first + 1:nb[b]
  
  num_zeros <- sum(y_obs[,ind_batch]==0)
  zero_rate[b] <- num_zeros/G/nb[b]
  num_dropouts <- sum(y_obs[,ind_batch]==0 & x_imputed[,ind_batch]>0)
  drop_rate[b] <- num_dropouts/num_zeros
  ind_first <- ind_first + nb[b]
}

print(zero_rate)
print(drop_rate)

save.image(paste0("Dropout_",proj,".RData"))
