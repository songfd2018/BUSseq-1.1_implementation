rm(list = ls())
# run BUSseq on the benchmarked dataset

library(scran)
library(scater)
library(ggplot2)
library(Rtsne)
library(RColorBrewer)
library(xtable)

setwd("RareCellType")
################
# load dataset #
################
load("../sc_mixology-master/data/9cellmix_qc.RData")

###################################################
# feasure selection to generate the observed data #
###################################################
Num_hvg <- 6000

proj <- "LUAD"
interested_type <- c("9-0-0","0-9-0","0-0-9")
# following the workflow of MNN
# as http://bioconductor.org/packages/3.7/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html

cell_num_sel <- c(10,8,6,5,4,3,2,1)
ver <- 4

for(n_sel in cell_num_sel){

## cell_mix1 dataset
## cell selection
celltype <- paste0(sce_SC1_qc$HCC827, "-", sce_SC1_qc$H1975, "-", sce_SC1_qc$H2228)
sce_SC1_filter <- sce_SC1_qc[,celltype %in% interested_type ]

# downsample the cell line (9-0-0) to regard it as the rare cell type
celltype_select <- celltype[ celltype %in% interested_type ]
first_celltype <- which(celltype_select=="9-0-0")

# sample the cells to be removed
cell_removed <- sample(first_celltype, length(first_celltype) - n_sel)
celltype_select <- celltype_select[-cell_removed]
sce_SC1_filter <- sce_SC1_filter[, -cell_removed]

## normlization
clusters <- quickCluster(sce_SC1_filter, min.mean=0.1)
sce_SC1_filter <- computeSumFactors(sce_SC1_filter, min.mean=0.1, clusters=clusters)
sce_SC1_filter <- normalize(sce_SC1_filter)

var.fit <- trendVar(sce_SC1_filter, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_SC1_filter, var.fit)

## select HVGs
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:Num_hvg], ]
hvg_mix1 <- rownames(hvg.out)

## cell_mix2 dataset
## cell selection
celltype <- paste0(sce_SC2_qc$HCC827, "-", sce_SC2_qc$H1975, "-", sce_SC2_qc$H2228)
sce_SC2_filter <- sce_SC2_qc[,celltype %in% interested_type ]

# downsample the cell line (9-0-0) to regard it as the rare cell type
celltype_select <- celltype[ celltype %in% interested_type ]
first_celltype <- which(celltype_select=="9-0-0")

# sample the cells to be removed
cell_removed <- sample(first_celltype, length(first_celltype) - n_sel)
celltype_select <- celltype_select[-cell_removed]
sce_SC2_filter <- sce_SC2_filter[, -cell_removed]

## normlization
clusters <- quickCluster(sce_SC2_filter, min.mean=0.1)
sce_SC2_filter <- computeSumFactors(sce_SC2_filter, min.mean=0.1, clusters=clusters)
sce_SC2_filter <- normalize(sce_SC2_filter)

var.fit <- trendVar(sce_SC2_filter, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_SC2_filter, var.fit)

## select HVGs
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:Num_hvg], ]
hvg_mix2 <- rownames(hvg.out)

## cell_mix3 dataset
## cell selection
celltype <- paste0(sce_SC3_qc$HCC827, "-", sce_SC3_qc$H1975, "-", sce_SC3_qc$H2228)
sce_SC3_filter <- sce_SC3_qc[,celltype %in% interested_type ]

# downsample the cell line (9-0-0) to regard it as the rare cell type
celltype_select <- celltype[ celltype %in% interested_type ]
first_celltype <- which(celltype_select=="9-0-0")

# sample the cells to be removed
cell_removed <- sample(first_celltype, length(first_celltype) - n_sel)
celltype_select <- celltype_select[-cell_removed]
sce_SC3_filter <- sce_SC3_filter[, -cell_removed]

## normlization
clusters <- quickCluster(sce_SC3_filter, min.mean=0.1)
sce_SC3_filter <- computeSumFactors(sce_SC3_filter, min.mean=0.1, clusters=clusters)
sce_SC3_filter <- normalize(sce_SC3_filter)

var.fit <- trendVar(sce_SC3_filter, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_SC3_filter, var.fit)

## select HVGs
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:Num_hvg], ]
hvg_mix3 <- rownames(hvg.out)

## cell_mix4 dataset
## cell selection
celltype <- paste0(sce_SC4_qc$HCC827, "-", sce_SC4_qc$H1975, "-", sce_SC4_qc$H2228)
sce_SC4_filter <- sce_SC4_qc[,celltype %in% interested_type ]

# downsample the cell line (9-0-0) to regard it as the rare cell type
celltype_select <- celltype[ celltype %in% interested_type ]
first_celltype <- which(celltype_select=="9-0-0")

# sample the cells to be removed
cell_removed <- sample(first_celltype, length(first_celltype) - n_sel)
celltype_select <- celltype_select[-cell_removed]
sce_SC4_filter <- sce_SC4_filter[, -cell_removed]

## normlization
clusters <- quickCluster(sce_SC4_filter, min.mean=0.1)
sce_SC4_filter <- computeSumFactors(sce_SC4_filter, min.mean=0.1, clusters=clusters)
sce_SC4_filter <- normalize(sce_SC4_filter)

var.fit <- trendVar(sce_SC4_filter, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_SC4_filter, var.fit)

## select HVGs
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:Num_hvg], ]
hvg_mix4 <- rownames(hvg.out)

# Combined HVG
common.hvg <- intersect(hvg_mix1, intersect(hvg_mix2, intersect(hvg_mix3, hvg_mix4)))

####################
# Combined dataset #
####################
# index of common HVG in mix1 dataset
match_index <- match(common.hvg, rownames(sce_SC1_filter))
CountMat_mix1 <- counts(sce_SC1_filter)[match_index,]

# index of common HVG in mix2 dataset
match_index <- match(common.hvg, rownames(sce_SC2_filter))
CountMat_mix2 <- counts(sce_SC2_filter)[match_index,]

# index of common HVG in mix3 dataset
match_index <- match(common.hvg, rownames(sce_SC3_filter))
CountMat_mix3 <- counts(sce_SC3_filter)[match_index,]

# index of common HVG in mix4 dataset
match_index <- match(common.hvg, rownames(sce_SC4_filter))
CountMat_mix4 <- counts(sce_SC4_filter)[match_index,]

nb <- c(ncol(CountMat_mix1), ncol(CountMat_mix2), ncol(CountMat_mix3), ncol(CountMat_mix4))

LUADmixCounts <- cbind(CountMat_mix1, CountMat_mix2, CountMat_mix3, CountMat_mix4)
colnames(LUADmixCounts) <- c(paste0("mix1_",colnames(CountMat_mix1)),
                             paste0("mix2_",colnames(CountMat_mix2)),
                             paste0("mix3_",colnames(CountMat_mix3)),
                             paste0("mix4_",colnames(CountMat_mix4)))

# metadata
Num_Cell1 <- c(sce_SC1_filter$HCC827, sce_SC2_filter$HCC827, sce_SC3_filter$HCC827, sce_SC4_filter$HCC827)
Num_Cell2 <- c(sce_SC1_filter$H1975, sce_SC2_filter$H1975, sce_SC3_filter$H1975, sce_SC4_filter$H1975)
Num_Cell3 <- c(sce_SC1_filter$H2228, sce_SC2_filter$H2228, sce_SC3_filter$H2228, sce_SC4_filter$H2228)

batch_label <- rep(c("Mix1","Mix2","Mix3","Mix4"), nb)

metadata <- data.frame(Batch = batch_label, 
                       Num_Cell1 = Num_Cell1,
                       Num_Cell2 = Num_Cell2,
                       Num_Cell3 = Num_Cell3)
rownames(metadata) <- colnames(LUADmixCounts)


#store count data matrix
if(!dir.exists("RawCountData")){
  dir.create("RawCountData")
}
write.table(LUADmixCounts,file=paste0("RawCountData/count_data_",proj,"_v",ver,".txt"),row.names = FALSE, col.names = FALSE)
write.table(metadata,file=paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

N <- ncol(LUADmixCounts)
G <- nrow(LUADmixCounts)
B <- length(nb)
write.table(c(N,G,B,nb),file=paste0("RawCountData/dim_",proj,"_v",ver,".txt"),row.names = FALSE, col.names = FALSE)
write.table(common.hvg,paste0("RawCountData/gene_list_",proj,"_v",ver,".txt"),row.names = FALSE, col.names = FALSE)

if(!dir.exists(paste0("./v",ver))){
  dir.create(paste0("./v",ver))
}

ver <- ver + 1
}