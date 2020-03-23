#Apply scanorama to the simulated study.
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
setwd("MainSimulation/Comparison/scanorama/")

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
scanorama_embedding <- read.table("scanorama_simulation_v1_integrated.txt")
scanorama_corrected <- read.table("scanorama_simulation_v1_corrected.txt")

clu_scanorama <- pam(scanorama_corrected, 7)
w_scanorama <- clu_scanorama$clustering

#######
# ARI #
#######
ARI_scanorama <- adjustedRandIndex(metadata$celltype, w_scanorama)

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
scanorama_dist <- dist(scanorama_embedding) 
all.dists.scanorama <- as.matrix(scanorama_dist)
tsne_scanorama_dist <- Rtsne(all.dists.scanorama, is_distance=TRUE, perplexity = 30)

scanorama_by_celltype<- "Image/tsne_simulation_scanorama_by_celltype.jpeg"
dat_frame <- data.frame(Var1 = tsne_scanorama_dist$Y[,1], 
                        Var2 = tsne_scanorama_dist$Y[,2], 
                        col = celltype_factor)

jpeg(scanorama_by_celltype,width = 900, height = 600, quality = 100)
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

scanorama_by_batch <- "Image/tsne_simulation_scanorama_by_batch.jpeg"
dat_frame <- data.frame(Var1 = tsne_scanorama_dist$Y[,1], 
                        Var2 = tsne_scanorama_dist$Y[,2], 
                        col = batch_factor)

jpeg(scanorama_by_batch,width = 900, height = 600, quality = 100)
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
pca.scanorama <- prcomp(scanorama_embedding, rank=2)

scanorama_PCA<- "Image/pca_simulation_scanorama.jpeg"
dat_frame <- data.frame(Var1 = pca.scanorama$x[,1], 
                        Var2 = pca.scanorama$x[,2], 
                        col = celltype_factor)
jpeg(scanorama_PCA,width = 900, height = 600, quality = 100)
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
save(ARI_scanorama,tsne_scanorama_dist,file = "scanorama_results.RData")