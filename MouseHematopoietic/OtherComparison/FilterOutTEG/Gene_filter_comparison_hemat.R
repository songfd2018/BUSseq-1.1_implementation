#Apply BUSseq to the hematopoietic study.
rm(list=ls())
library(mclust) # For ARI
library(xtable)
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

proj <- "hemat"
ver_ori <- 1
ver_cur <- 3

set.seed(123)
###################
# Load hemat Data #
###################
# Working directory
setwd("OtherComparison/FilterOutTEG")

# Loading hemat count data
y_ori <- read.table(paste0("../../RawCountData/count_data_",proj,"_v",ver_ori,".txt"))
y_cur <- read.table(paste0("RawCountData/count_data_",proj,"_v",ver_cur,".txt"))

# Load dimension
dim_ori <- read.table(paste0("../../RawCountData/dim_",proj,"_v",ver_ori,".txt"))
dim_ori <- unlist(dim_ori)
N <- dim_ori[1]
G_ori <- dim_ori[2]
B <- dim_ori[3]
nb <- dim_ori[3+1:B]

dim_cur <- read.table(paste0("RawCountData/dim_",proj,"_v",ver_cur,".txt"))
dim_cur <- unlist(dim_cur)
G_cur <- dim_cur[2]

# Load metadata
metadata <- read.table(paste0("../../RawCountData/metadata_",proj,"_v",ver_ori,".txt"))

# Load gene_list
gene_list_ori <- unlist(read.table(paste0("../../RawCountData/gene_list_",proj,"_v",ver_ori,".txt"),stringsAsFactors = F))
gene_list_cur <- unlist(read.table(paste0("RawCountData/gene_list_",proj,"_v",ver_cur,".txt"),stringsAsFactors = F))


##################################
# Compare the cell type labeling #
##################################
K <- 6

# load w_est
w_ori <- read.table(paste0("../../Inference_K",K,"/w_est.txt"))
w_ori <- unlist(w_ori)

w_cur <- read.table(paste0("Inference_K",K,"/w_est.txt"))
w_cur <- unlist(w_cur)

table(metadata$CellType, w_ori)[c(4,6,3,7,1,2,5),c(6,4,1,2,3,5)]
table(metadata$CellType, w_cur)[c(4,6,3,7,1,2,5),c(6,4,1,2,5,3)]
Com_oc <- table(w_cur, w_ori)[c(6,4,1,2,5,3),c(6,4,1,2,3,5)]
colnames(Com_oc) <- paste0("Cluster",1:K)
rownames(Com_oc) <- paste0("Cluster",1:K)
xtable(Com_oc)

###########################
# Compare intrinsic genes #
###########################
# load intrinsic genes
D_ori <- unlist(read.table(paste0("../../Inference_K",K,"/IG_est.txt")))
intri_ori <- gene_list_ori[D_ori==1]

D_cur <- unlist(read.table(paste0("Inference_K",K,"/IG_est.txt")))
intri_cur <- gene_list_cur[D_cur==1]

shared_genes <- intersect(gene_list_ori, gene_list_cur)

intri_share_ori <- intersect(shared_genes, intri_ori) 
intri_share_cur <- intersect(shared_genes, intri_cur)

length(intri_share_ori)
length(intri_share_cur)
length(intersect(intri_share_ori,intri_share_cur))

save.image("Comparison_HVG_HEG_hemat.RData")