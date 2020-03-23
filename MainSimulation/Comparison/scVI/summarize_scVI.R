# Apply scVI to the simulated study.
rm(list=ls())
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
setwd("MainSimulation/Comparison/scVI/")

###################
# Load simulation Data #
###################
# load the dimentsion information
dim <- unlist(read.table("../../RawCountData/dim_simulation_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[1:B + 3]

# Loading the true cell type indicators
metadata <- read.table("../../RawCountData/metadata_simulation_v1.txt")

##############################
# load the inference by scVI #
##############################
w_scVI <- unlist(read.table("scVI_simulation_v1_clusters.txt"))
scVI_corrected <- read.table("scVI_simulation_v1_latent.txt")

#######
# ARI #
#######
ARI_scVI <- adjustedRandIndex(metadata$celltype, w_scVI)

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
set.seed(123)
scVI_dist <- dist(scVI_corrected) 
all.dists.scVI <- as.matrix(scVI_dist)
tsne_scVI_dist <- Rtsne(all.dists.scVI, is_distance=TRUE, perplexity = 30)

scVI_by_celltype<- "Image/tsne_simulation_scVI_by_celltype.jpeg"
dat_frame <- data.frame(Var1 = tsne_scVI_dist$Y[,1], 
                        Var2 = tsne_scVI_dist$Y[,2], 
                        col = celltype_factor)

jpeg(scVI_by_celltype,width = 900, height = 600, quality = 100)
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

scVI_by_batch <- "Image/tsne_simulation_scVI_by_batch.jpeg"
dat_frame <- data.frame(Var1 = tsne_scVI_dist$Y[,1], 
                        Var2 = tsne_scVI_dist$Y[,2], 
                        col = batch_factor)

jpeg(scVI_by_batch,width = 900, height = 600, quality = 100)
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

##################
# Draw PCA plots #
##################
pca.scVI <- prcomp(scVI_corrected, rank=2)

scVI_PCA<- "Image/pca_simulation_scVI.jpeg"
dat_frame <- data.frame(Var1 = pca.scVI$x[,1], 
                        Var2 = pca.scVI$x[,2], 
                        col = celltype_factor)
jpeg(scVI_PCA,width = 900, height = 600, quality = 100)
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
save(ARI_scVI,tsne_scVI_dist,file = "scVI_results.RData")