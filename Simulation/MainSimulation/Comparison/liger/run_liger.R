# library(devtools)
# install_github('MacoskoLab/liger')
rm(list=ls())
library(liger)
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
setwd("MainSimulation/Comparison/liger/")

###################
# Load simulation Data #
###################
# Loading simulation count data
countdata <- read.table("../../RawCountData/count_data_simulation_v1.txt") 

# Load dimension
dim <- unlist(read.table("../../RawCountData/dim_simulation_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_simulation_v1.txt")

# Load gene_list
gene_list <- paste0("gene_",1:G)

rownames(countdata) <- gene_list
colnames(countdata) <- rownames(metadata)
simulation.list <- list(Batch1 = countdata[,1:nb[1]], 
                        Batch2 = countdata[,1:nb[2]+nb[1]], 
                        Batch3 = countdata[,1:nb[3]+sum(nb[1:2])], 
                        Batch4 = countdata[,1:nb[4]+sum(nb[1:3])])

################################################################################################################
# Conduct analysis referring to https://github.com/MacoskoLab/liger/blob/master/vignettes/walkthrough_pbmc.pdf #
################################################################################################################
# Create liger object
liger_simulation <- createLiger(simulation.list, remove.missing = FALSE)

# Normalization
liger_simulation <- normalize(liger_simulation)
liger_simulation@var.genes <- gene_list
liger_simulation <- scaleNotCenter(liger_simulation)

# running suggestK on multiple cores can greatly decrease the runtime
# looking for the section of the plot where this metric stops increasing as sharply
k.suggest <- suggestK(liger_simulation, num.cores = 5, gen.new = T, return.data = T, plot.log2 = F)


# Take the lowest objective of three factorizations with different initializations
# Multiple restarts are recommended for initial analyses since iNMF is non-deterministic
liger_simulation <- optimizeALS(liger_simulation, k=30, thresh = 5e-5, nrep = 3)

liger_simulation <- runTSNE(liger_simulation, use.raw = T)
p1 <- plotByDatasetAndCluster(liger_simulation, return.plots = T)
# Plot by dataset
print(p1[[1]])


##########################
# Currently working here #
##########################

# To better integrate the datasets, we perform a quantile alignment step. This process first identifies similarly
# loading cells across datasets by building a similarity graph based on shared factor neighborhoods.
liger_simulation <- quantileAlignSNF(liger_simulation, resolution = 0.4, small.clust.thresh = 20)
w_liger <- liger_simulation@alignment.clusters

#######
# ARI #
#######
ARI_liger <- adjustedRandIndex(metadata$celltype,w_liger)

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
loading_liger <- liger_simulation@H.norm
liger_dist <- dist(loading_liger)
all.dists.liger <- as.matrix(liger_dist)
tsne_liger_dist <- Rtsne(all.dists.liger, is_distance=TRUE, perplexity = 30)

liger_by_celltype<- "Image/tsne_simulation_liger_by_celltype.jpeg"
dat_frame <- data.frame(Var1 = tsne_liger_dist$Y[,1], 
                        Var2 = tsne_liger_dist$Y[,2], 
                        col = celltype_factor)

jpeg(liger_by_celltype,width = 900, height = 600, quality = 100)
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

liger_by_batch <- "Image/tsne_simulation_liger_by_batch.jpeg"
dat_frame <- data.frame(Var1 = tsne_liger_dist$Y[,1], 
                        Var2 = tsne_liger_dist$Y[,2], 
                        col = batch_factor)

jpeg(liger_by_batch,width = 900, height = 600, quality = 100)
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
pca.liger <- prcomp(loading_liger, rank=2)

liger_PCA<- "Image/pca_simulation_liger.jpeg"
dat_frame <- data.frame(Var1 = pca.liger$x[,1], 
                        Var2 = pca.liger$x[,2], 
                        col = celltype_factor)
jpeg(liger_PCA,width = 900, height = 600, quality = 100)
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
save(ARI_liger,tsne_liger_dist,file = "liger_results.RData")