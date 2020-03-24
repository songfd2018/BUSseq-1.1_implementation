# Apply Seurat 3 to the pancreatic study.
# Referring to https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html
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
setwd("BatchEffectCorrectionComparison/Seurat/")

########################
# Load pancreatic Data #
########################

# Loading hematopoietic count data
countdata <- read.table("../../RawCountData/count_data_pancreas_v1.txt") 

dim <- unlist(read.table("../../RawCountData/dim_pancreas_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_pancreas_v1.txt")

# Load gene_list
gene_list <- unlist(read.table("../../RawCountData/gene_list_pancreas_v1.txt",stringsAsFactors = F))

rownames(countdata) <- gene_list
colnames(countdata) <- rownames(metadata)

##########################################
# Apply Seurat to the Hematopoietic Data #
##########################################
# Setting the Seurat Object
pancreas <- CreateSeuratObject(countdata, meta.data = metadata)
pancreas.list <- SplitObject(pancreas, split.by = "Study")


for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = G, 
                                             verbose = FALSE)
}

reference.list <- pancreas.list[c("GSE81076", "GSE85241","GSE86473","E-MTAB-5061")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = G)

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

DefaultAssay(pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)

pancreas.integrated <- FindNeighbors(pancreas.integrated, reduction = "pca", dims = 1:30, k.param = 20)
pancreas.integrated<-FindClusters(pancreas.integrated,resolution = 1.5)

w_Seurat <- pancreas.integrated$seurat_clusters
reorder <- match(rownames(metadata),names(w_Seurat))
w_Seurat<- w_Seurat[reorder]

Seurat_PCA <- pancreas.integrated@reductions$pca@cell.embeddings[reorder,]

Seurat_corrected <- NULL
for(g in 1:G){
  Seurat_corrected <- rbind(Seurat_corrected, pancreas.integrated@assays$integrated[g,])
}

#######
# ARI #
#######
ARI_Seurat <- adjustedRandIndex(metadata$CellType,w_Seurat)

kclust_Seurat <- pam(pancreas.integrated@reductions$pca@cell.embeddings, k = 7)
w_Seurat_kmeans <- kclust_Seurat$clustering
ARI_Seurat_kmeans <- adjustedRandIndex(w_Seurat_kmeans,metadata$CellType)

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
Seurat_dist <- dist(Seurat_PCA)
all.dists.Seurat <- as.matrix(Seurat_dist)
tsne_Seurat_dist <- Rtsne(all.dists.Seurat, is_distance=TRUE, perplexity = 30)

Seurat_by_celltype<- "Image/tsne_pancreas_Seurat_by_celltype.jpg"
dat_frame <- data.frame(Var1 = tsne_Seurat_dist$Y[,1], 
                        Var2 = tsne_Seurat_dist$Y[,2], 
                        col = celltype_factor)

jpeg(Seurat_by_celltype,width = 800, height = 600, quality = 100)
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

Seurat_by_batch <- "Image/tsne_pancreas_Seurat_by_batch.jpg"
dat_frame <- data.frame(Var1 = tsne_Seurat_dist$Y[,1], 
                        Var2 = tsne_Seurat_dist$Y[,2], 
                        col = batch_factor)

jpeg(Seurat_by_batch,width = 800, height = 600, quality = 100)
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
Seu_PCA<- "Image/pca_pancreas_Seurat.jpeg"
dat_frame <- data.frame(Var1 = pca.Seurat$x[,1], 
                        Var2 = pca.Seurat$x[,2], 
                        col = celltype_factor)

jpeg(Seu_PCA,width = 800, height = 600, quality = 100)
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
save(ARI_Seurat, ARI_Seurat_kmeans, tsne_Seurat_dist, Seurat_PCA, file = "Seurat_results.RData")