#Apply MNN to the simualtion study.
rm(list=ls())
library(scran)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)

# coloring
library(scales)
require(WGCNA)
library(RColorBrewer)


set.seed(12345)
# Working directory
setwd("MainSimulation/Comparison/MNN/")

########################
# Load simulation Data #
########################
# Loading the file name list of all simulation count data
countdata <- read.table("../../RawCountData/count_data_simulation_v1.txt") 

dim <- unlist(read.table("../../RawCountData/dim_simulation_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_simulation_v1.txt")
gene_list <- paste0("gene_",1:G)

colnames(countdata) <- rownames(metadata)
rownames(countdata) <- gene_list

simulationCounts <- list()
index <- 0
for(b in 1:B){
  simulationCounts[[b]] <- as.matrix(countdata[,index + 1:nb[b]])
  index <- index + nb[b]
}

##################################
# Apply MNN to the simulation Data #
##################################
data_MNN <- simulationCounts
# Referring to https://github.com/MarioniLab/MNN2017/blob/master/simulation
# Referring to http://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/batch.html
# Normalization
data_MNN_normalized<-list()

for(b in 1:B){
  high.abF <- scater::calcAverage(data_MNN[[b]]) > 1
  clustF <- quickCluster(data_MNN[[b]], min.size=10 , method="igraph", subset.row=high.abF)
  sizeF <- computeSumFactors(data_MNN[[b]], sizes=seq(11, 81, 5), cluster=clustF, subset.row=high.abF)
  data_MNN_normalized[[b]] <- t(t(data_MNN[[b]])/sizeF)
}

# Rescaling the first dataset to match the coverage of the second.
ave<-list()
for(b in 1:B){
  ave[[b]] <- rowMeans(data_MNN_normalized[[b]])
  ave[[b]] <- max(1e-6,ave[[b]])
  if(b>1){
    data_MNN_normalized[[b]] <- data_MNN_normalized[[b]] * median(ave[[1]]/ave[[b]])
  }
}

# MNN batch correction
Xmnn <- mnnCorrect(data_MNN_normalized[[1]],
                   data_MNN_normalized[[2]],
                   data_MNN_normalized[[3]],
                   data_MNN_normalized[[4]],
                   svd.dim=0,
                   cos.norm.in=TRUE, cos.norm.out=TRUE,
                   var.adj=TRUE, 
                   k=20, sigma=0.1)

# combine corrected matrices together
corrected.df <- do.call(cbind.data.frame, Xmnn$corrected)
corrected.mat <- as.matrix(t(corrected.df))

##############
# Clustering #
##############

omat <- do.call(cbind, data_MNN_normalized)
sce <- SingleCellExperiment(list(logcounts=log1p(omat)))
reducedDim(sce, "Corrected") <- corrected.mat
sce$Batch <- rep(paste0("Batch",1:B),nb)

start_time <- Sys.time()
snn.gr <- buildSNNGraph(sce, use.dimred="Corrected")
clusters <- igraph::cluster_walktrap(snn.gr)
end_time <- Sys.time()
time_consumption <- end_time - start_time

table(clusters$membership, sce$Batch)

# The estimated cell type indicators by MNN
w_MNN <- factor(clusters$membership)

#######
# ARI #
#######
ARI_MNN <- adjustedRandIndex(metadata$celltype,w_MNN)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}

#set cell type colorings
# batch color bar
color_by_batch<-c("#EB4334","#FBBD06","#35AA53","#4586F3")
batch_factor <- factor(paste0("Batch ",rep(1:B,nb)))

# cell type color bar
color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978","#8097D3","#9C7ACE")
celltype_factor <- factor(paste0("Celltype ",metadata$celltype))

#################################################
# Draw t-SNE plots by batch and true cell types #
#################################################
MNN_dist <- dist(corrected.mat)
set.seed(123)
all.dists.MNN <- as.matrix(MNN_dist)
tsne_MNN_dist <- Rtsne(all.dists.MNN, is_distance=TRUE, perplexity = 30)

MNN_by_celltype<- "Image/tsne_simulation_MNN_by_celltype.jpeg"
dat_frame <- data.frame(Var1 = tsne_MNN_dist$Y[,1], 
                        Var2 = tsne_MNN_dist$Y[,2], 
                        col = celltype_factor)

jpeg(MNN_by_celltype,width = 900, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_celltype,0.6)) +
  xlab("tSNE 1") + ylab("tSNE 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 44), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 44),#, angle = 45))
        axis.title=element_text(size=48,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.position = "none")
dev.off()

MNN_by_batch <- "Image/tsne_simulation_MNN_by_batch.jpeg"
dat_frame <- data.frame(Var1 = tsne_MNN_dist$Y[,1], 
                        Var2 = tsne_MNN_dist$Y[,2], 
                        col = batch_factor)

jpeg(MNN_by_batch,width = 900, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_batch,0.6)) +
  xlab("tSNE 1") + ylab("tSNE 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.position = "none")
dev.off()

##########################################
# Draw the PCA plot on common cell types #
##########################################
pca.MNN <- prcomp(corrected.mat, rank=2)

MNN_PCA<- "Image/pca_simulation_MNN.jpeg"
dat_frame <- data.frame(Var1 = pca.MNN$x[,1], 
                        Var2 = pca.MNN$x[,2], 
                        col = celltype_factor)
jpeg(MNN_PCA,width = 900, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_celltype,0.6)) +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.position = "none")
dev.off()

# Store the workspace
save(ARI_MNN,tsne_MNN_dist,file = "MNN_results.RData")