# Apply Seurat 3 to the simulated study.
# Referring to https://satijalab.org/seurat/v3.0/simulation_integration_label_transfer.html
rm(list=ls())
library(Seurat)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)
library(cowplot)

# coloring
library(scales)
require(WGCNA)
library(RColorBrewer)

set.seed(123)
# Working directory
setwd("MainSimulation/Comparison/Seurat/")

########################
# Load simulated Data #
########################

# Loading hematopoietic count data
countdata <- read.table("../../RawCountData/count_data_simulation_v1.txt") 

dim <- unlist(read.table("../../RawCountData/dim_simulation_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_simulation_v1.txt")

# Load gene_list
gene_list <- paste0("gene-",1:G)

rownames(countdata) <- gene_list
colnames(countdata) <- rownames(metadata)

##########################################
# Apply Seurat to the Simulated Data #
##########################################
# Setting the Seurat Object
simulation <- CreateSeuratObject(countdata, meta.data = metadata)
simulation.list <- SplitObject(simulation, split.by = "batch")


for (i in 1:length(simulation.list)) {
  simulation.list[[i]] <- NormalizeData(simulation.list[[i]], verbose = FALSE)
  simulation.list[[i]] <- FindVariableFeatures(simulation.list[[i]], selection.method = "vst", nfeatures = G, 
                                             verbose = FALSE)
}

reference.list <- simulation.list[c("Batch_1", "Batch_2","Batch_3","Batch_4")]
simulation.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = G)

simulation.integrated <- IntegrateData(anchorset = simulation.anchors, dims = 1:30)

DefaultAssay(simulation.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
simulation.integrated <- ScaleData(simulation.integrated, verbose = FALSE)
simulation.integrated <- RunPCA(simulation.integrated, npcs = 30, verbose = FALSE)

# simulation.integrated <- RunUMAP(simulation.integrated, reduction = "pca", dims = 1:30)
# 
# p1 <- DimPlot(simulation.integrated, reduction = "umap", group.by = "Study")
# p2 <- DimPlot(simulation.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
#   NoLegend()
# plot_grid(p1, p2)

simulation.integrated <- FindNeighbors(simulation.integrated, reduction = "pca", dims = 1:30, k.param = 20)
simulation.integrated<-FindClusters(simulation.integrated,resolution = 1.5)

# simulation.integrated <- RunTSNE(object = simulation.integrated)
# DimPlot(object = simulation.integrated, reduction = "tsne")

w_Seurat <- simulation.integrated$seurat_clusters
reorder <- match(rownames(metadata),names(w_Seurat))
w_Seurat<- w_Seurat[reorder]

Seurat_PCA <- simulation.integrated@reductions$pca@cell.embeddings[reorder,]

Seurat_corrected <- NULL
for(g in 1:G){
  Seurat_corrected <- rbind(Seurat_corrected, simulation.integrated@assays$integrated[g,])
}

#######
# ARI #
#######
ARI_Seurat <- adjustedRandIndex(metadata$celltype,w_Seurat)

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
Seurat_dist <- dist(Seurat_PCA)
all.dists.Seurat <- as.matrix(Seurat_dist)
tsne_Seurat_dist <- Rtsne(all.dists.Seurat, is_distance=TRUE, perplexity = 30)

Seurat_by_celltype<- "Image/tsne_simulation_Seurat_by_celltype.jpeg"
dat_frame <- data.frame(Var1 = tsne_Seurat_dist$Y[,1], 
                        Var2 = tsne_Seurat_dist$Y[,2], 
                        col = celltype_factor)

jpeg(Seurat_by_celltype,width = 900, height = 600, quality = 100)
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

Seurat_by_batch <- "Image/tsne_simulation_Seurat_by_batch.jpeg"
dat_frame <- data.frame(Var1 = tsne_Seurat_dist$Y[,1], 
                        Var2 = tsne_Seurat_dist$Y[,2], 
                        col = batch_factor)

jpeg(Seurat_by_batch,width = 900, height = 600, quality = 100)
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
pca.Seurat <- prcomp(Seurat_PCA, rank=2)
Seu_PCA<- "Image/pca_simulation_Seurat.jpeg"
dat_frame <- data.frame(Var1 = pca.Seurat$x[,1], 
                        Var2 = pca.Seurat$x[,2], 
                        col = celltype_factor)
jpeg(Seu_PCA,width = 900, height = 600, quality = 100)
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
save(ARI_Seurat,tsne_Seurat_dist,file = "Seurat_results.RData")

