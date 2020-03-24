# Apply ZINB-WaVE to the pancreas study.
# Referring to https://github.com/drisso/zinb_analysis/blob/master/real_data/allen_covariates_1000.Rmd
rm(list=ls())
library(zinbwave)
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
setwd("BatchEffectCorrectionComparison/ZINBWaVE/")

########################
# Load Pancreatic Data #
########################
# Loading the file name list of all pancreas count data
countdata <- read.table("../../RawCountData/count_data_pancreas_v1.txt") 

dim <- unlist(read.table("../../RawCountData/dim_pancreas_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_pancreas_v1.txt")
gene_list <- unlist(read.table("../../RawCountData/gene_list_pancreas_v1.txt",stringsAsFactors = FALSE))

colnames(countdata) <- rownames(metadata)
rownames(countdata) <- gene_list


##########################################
# Apply ZINB-WaVE to the Pancreatic Data #
##########################################
data_ZINBW <- as.matrix(countdata)

# Factorizing the batch indicators
batch_ind <- factor(rep(1:B,nb))

# Performing ZINB-WaVE
zinb_batch <- zinbFit(data_ZINBW, K = 10, X=model.matrix(~batch_ind), epsilon=1e3)

# Clustering
data_ZINBW <- SummarizedExperiment(assays = data_ZINBW)
assayNames(data_ZINBW)[1] <- "counts"
merged_zinb <- zinbwave(data_ZINBW, fitted_model = zinb_batch, K = 10, epsilon=1000)

library(Seurat)

seu <- as.Seurat(x = merged_zinb, counts = "counts", data = "counts")

seu
seu <- FindNeighbors(seu, reduction = "zinbwave",
                     dims = 1:10 #this should match K
)
seu <- FindClusters(object = seu)

w_ZINBWaVE <- seu$seurat_clusters

#######
# ARI #
#######
ARI_ZINBWaVE <- adjustedRandIndex(metadata$CellType,w_ZINBWaVE)

kclust_ZINBWaVE <- pam(zinb_batch@W, k = 7)
w_ZINBWaVE_kmeans <- kclust_ZINBWaVE$clustering
ARI_kmeans_kmeans <- adjustedRandIndex(w_ZINBWaVE_kmeans,metadata$CellType)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}

#####set cell type colorings
color_by_celltype <- c(brewer.pal(6,"Set3"),"black")
color_by_celltype[c(1,6)] <- color_by_celltype[c(6,1)]
celltype_factor <- factor(metadata$CellType, levels = c("Alpha", "Beta", "Gamma", "Delta", "Acinar", "Ductal", "other"))

#####set batch colorings
color_by_batch <- c("#EB4334","#FBBD06","#35AA53","#4586F3")
batch_factor <- factor(rep(1:B,nb))

#################################################
# Draw t-SNE plots by batch and true cell types #
#################################################
set.seed(123)
ZINBWaVE_dist <- dist(zinb_batch@W)
all.dists.ZINBW <- as.matrix(ZINBWaVE_dist)
tsne_ZINBW_dist <- Rtsne(all.dists.ZINBW, is_distance=TRUE, perplexity = 30)

ZINBW_by_celltype<- "Image/tsne_pancreas_ZINBW_by_celltype.jpg"
dat_frame <- data.frame(Var1 = tsne_ZINBW_dist$Y[,1], 
                        Var2 = tsne_ZINBW_dist$Y[,2], 
                        col = celltype_factor)

jpeg(ZINBW_by_celltype,width = 800, height = 600, quality = 100)
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

ZINBW_by_batch <- "Image/tsne_pancreas_ZINBW_by_batch.jpg"
dat_frame <- data.frame(Var1 = tsne_ZINBW_dist$Y[,1], 
                        Var2 = tsne_ZINBW_dist$Y[,2], 
                        col = batch_factor)

jpeg(ZINBW_by_batch,width = 800, height = 600, quality = 100)
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
pca.ZINBW <- prcomp(zinb_batch@W, rank=2)
ZINBW_PCA<- "./Image/pca_pancreas_ZINBWaVE.jpeg"
dat_frame <- data.frame(Var1 = pca.ZINBW$x[,1], 
                        Var2 = pca.ZINBW$x[,2], 
                        col = celltype_factor)

jpeg(ZINBW_PCA,width = 800, height = 600, quality = 100)
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
save(ARI_ZINBWaVE, ARI_kmeans_kmeans, tsne_ZINBW_dist, zinb_batch, file = "ZINBWaVE_results.RData")