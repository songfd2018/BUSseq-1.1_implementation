# Apply Seurat 3 to the hematopoietic study.
rm(list=ls())
library(Seurat)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
set.seed(123)


# Working directory
setwd("BatchEffectCorrectionComparison/Seurat/")

###################
# Load hemat Data #
###################
# Loading hematopoietic count data
countdata <- read.table("../../RawCountData/count_data_hemat_v1.txt") 

dim <- unlist(read.table("../../RawCountData/dim_hemat_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_hemat_v1.txt")
rownames(metadata) <- paste0("Batch_",rep(c(1,2),nb),"_Cell_",c(1:nb[1],1:nb[2]))
gene_list <- unlist(read.table("../../RawCountData/gene_list_hemat_v1.txt",stringsAsFactors = FALSE))

colnames(countdata) <- paste0("Batch_",rep(c(1,2),nb),"_Cell_",c(1:nb[1],1:nb[2]))
rownames(countdata) <- gene_list

##########################################
# Apply Seurat to the Hematopoietic Data #
##########################################
# Referrring to https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html
# Setting the Seurat Object 
hemat <- CreateSeuratObject(countdata, meta.data = metadata)
hemat.list <- SplitObject(hemat, split.by = "Study")

for (i in 1:length(hemat.list)) {
  hemat.list[[i]] <- NormalizeData(hemat.list[[i]], verbose = FALSE)
  hemat.list[[i]] <- FindVariableFeatures(hemat.list[[i]], selection.method = "vst", nfeatures = 3470, 
                                             verbose = FALSE)
}

reference.list <- hemat.list[c("GSE72857", "GSE81682")]
hemat.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = 3470)

hemat.integrated <- IntegrateData(anchorset = hemat.anchors, dims = 1:30)

DefaultAssay(hemat.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
hemat.integrated <- ScaleData(hemat.integrated, verbose = FALSE)
hemat.integrated <- RunPCA(hemat.integrated, npcs = 30, verbose = FALSE)

hemat.integrated <- FindNeighbors(hemat.integrated, reduction = "pca", dims = 1:30, k.param = 20)
hemat.integrated<-FindClusters(hemat.integrated,resolution = 1.5)

w_Seurat <- hemat.integrated$seurat_clusters

Seurat_PCA <- hemat.integrated@reductions$pca@cell.embeddings
Seurat_corrected <- NULL
for(g in 1:G){
  Seurat_corrected <- rbind(Seurat_corrected, hemat.integrated@assays$integrated[g,])
}

#######
# ARI #
#######
ARI_Seurat <- adjustedRandIndex(metadata$CellType,w_Seurat)

kclust_Seurat <- pam(hemat.integrated@reductions$pca@cell.embeddings, k = 7)
w_Seurat_kmeans <- kclust_Seurat$clustering
ARI_Seurat_kmeans <- adjustedRandIndex(w_Seurat_kmeans,metadata$CellType)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}

#####set cell type colorings
# batch color bar
color_by_batch<-c("#EB4334","#4586F3")
batch_factor <- factor(rep(1:B,nb))

# cell type color bar
color_by_celltype<-c("#1A9850", # LTHSC
                     "#66C2A5", # MPP
                     "#4393C3", # LMPP 
                     "#FEB24C", # CMP
                     "#D7301F", # MEP
                     "#FFFF99", # GMP
                     "black" )# other

celltype_factor <- factor(metadata$CellType, levels = c("LTHSC", "MPP", "LMPP", "CMP", "MEP", "GMP", "other"))

#################################################
# Draw t-SNE plots by batch and true cell types #
#################################################
set.seed(123)
Seurat_dist <- dist(Seurat_PCA)
all.dists.Seurat <- as.matrix(Seurat_dist)
tsne_Seurat_dist <- Rtsne(all.dists.Seurat, is_distance=TRUE, perplexity = 30)

Seurat_by_celltype<- "Image/tsne_hemat_Seurat_by_celltype.jpg"
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

Seurat_by_batch <- "Image/tsne_hemat_Seurat_by_batch.jpg"
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

Seu_PCA<- "Image/pca_hemat_Seurat.jpeg"
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