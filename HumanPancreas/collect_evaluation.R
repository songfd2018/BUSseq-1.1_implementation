rm(list=ls())
library(mclust)
library(ggplot2)
library(reshape2)
library(Rtsne)
library(cluster)

###################
# Load results #
###################
# Working directory
setwd("/your/working/directory/BUSseq_implementation-1.0/MouseHematopoietic/")

load("BUSseq/BUSseq_results.RData")
comparison_list <- c("liger", "MNN", "scanorama", "scVI", "Seurat", "ZINBWaVE")
num_comparison <- length(comparison_list)

for(m in 1:num_comparison){
  load(paste0("Comparison/",comparison_list[m],"/",comparison_list[m],"_results.RData"))
}

# load the dimentsion information
dim <- unlist(read.table("RawCountData/dim_hemat_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[1:B + 3]

# Load metadata
metadata <- read.table("RawCountData/metadata_hemat_v1.txt")

# Load gene_list
gene_list <- unlist(read.table("RawCountData/gene_list_hemat_v1.txt",stringsAsFactors = F))

# Load raw count data 
data_raw <- as.matrix(read.table("RawCountData/count_data_hemat_v1.txt"))

#############################################################
# Evaluation of the Performance of All the Methods          #
# according to ARI, silhouette coefficients and t-SNE plots #
#############################################################
if(!dir.exists("Results")){
  dir.create("Results")
}

if(!dir.exists("Image")){
  dir.create("Image")
}

#######
# ARI #
#######
ARI_values <- rep(NA, num_comparison + 1)
names(ARI_values) <- c("BUSseq", comparison_list)


# Calculating the ARI
ARI_values[1] <- ARI_BUSseq
ARI_values[2] <- ARI_liger
ARI_values[3] <- ARI_MNN
ARI_values[4] <- ARI_scanorama
ARI_values[5] <- ARI_scVI
ARI_values[6] <- ARI_Seurat
ARI_values[7] <- ARI_ZINBWaVE

write.table(ARI_values, "Results/ARI_values.txt",col.names = F)

# Calculating the ARI by robust k-means clustering
ARI_kmeans <- rep(NA, num_comparison + 1)
names(ARI_kmeans) <- c("BUSseq", comparison_list)

ARI_kmeans[1] <- ARI_BUSseq_kmeans
ARI_kmeans[2] <- ARI_liger_kmeans
ARI_kmeans[3] <- ARI_MNN_kmeans
ARI_kmeans[4] <- ARI_scanorama_kmeans
ARI_kmeans[5] <- ARI_scVI_kmeans
ARI_kmeans[6] <- ARI_Seurat_kmeans
ARI_kmeans[7] <- ARI_ZINBWaVE_kmeans

write.table(ARI_kmeans, "Results/ARI_kmeans.txt",col.names = F)

###############################################################################
# Silhouette coefficient of tSNE by true cell type labels of each method #
###############################################################################
# Storing the Silhouette coefficients
sil_tSNE_true <- matrix(NA,N,num_comparison+1)
colnames(sil_tSNE_true) <- names(ARI_values)

BUSseq_dist_tSNE <- dist(tsne_BUSseq_dist$Y)
sil_tSNE_true[,1] <- silhouette(as.integer(metadata$CellType), dist = BUSseq_dist_tSNE)[,3]

liger_dist_tSNE <- dist(tsne_liger_dist$Y)
sil_tSNE_true[,2] <- silhouette(as.integer(metadata$CellType), dist = liger_dist_tSNE)[,3]

MNN_dist_tSNE <- dist(tsne_MNN_dist$Y)
sil_tSNE_true[,3] <- silhouette(as.integer(metadata$CellType), dist = MNN_dist_tSNE)[,3]

scanorama_dist_tSNE <- dist(tsne_scanorama_dist$Y)
sil_tSNE_true[,4] <- silhouette(as.integer(metadata$CellType), dist = scanorama_dist_tSNE)[,3]

scVI_dist_tSNE <- dist(tsne_scVI_dist$Y)
sil_tSNE_true[,5] <- silhouette(as.integer(metadata$CellType), dist = scVI_dist_tSNE)[,3]

Seurat_dist_tSNE <- dist(tsne_Seurat_dist$Y)
sil_tSNE_true[,6] <- silhouette(as.integer(metadata$CellType), dist = Seurat_dist_tSNE)[,3]

ZINBWaVE_dist_tSNE <- dist(tsne_ZINBW_dist$Y)
sil_tSNE_true[,7] <- silhouette(as.integer(metadata$CellType), dist = ZINBWaVE_dist_tSNE)[,3]

write.table(sil_tSNE_true,"./Results/Silhouette_coefs.txt",row.names = F)

# Drawing the boxplot plot
colnames(sil_tSNE_true) <- c("BUSseq", "LIGER", "MNN", "Seurat", "Scanorama", "scVI", "ZINB-WaVE")
sil_com_tSNE_true<-melt(sil_tSNE_true)
colnames(sil_com_tSNE_true) <- c("No_Cell","Method","Silhouette_Coefficient")

jpeg("Image/pancreas_boxplot_of_Silhouette_coefs_tSNE.jpg",width = 900, height = 600,quality = 100)
p <- ggplot(sil_com_tSNE_true, aes(x=Method, y=Silhouette_Coefficient,fill=Method)) +
  geom_boxplot() + xlab(NULL) + ylab(NULL)+theme_bw() +
  theme(text = element_text(size=48),
        axis.text.x =element_text(face = "bold", #color = "#993333", 
                                  size = 36,angle = 90), 
        axis.text.y = element_text(face = "bold", #color = "#993333", 
                                   size = 36),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

# Storing the Silhouette coefficients
sil_pca_true <- matrix(NA,N,7)
pc_rank <- 10
colnames(sil_pca_true) <- c("BUSseq", "LIGER", "MNN", "Scanorama", "scVI", "Seurat", "ZINBWaVE")

# 10-dimension PCs
PCs_BUSseq <- prcomp(t(log1p(x_intrinsic)), rank=pc_rank)$x
BUSseq_dist_pca <- dist(PCs_BUSseq)
sil_pca_true[,1] <- silhouette(as.integer(metadata$CellType), dist = BUSseq_dist_pca)[,3]

# load the workspace of liger
PCs_liger <- prcomp(loading_liger, rank=pc_rank)$x
liger_dist_pca <- dist(PCs_liger)
sil_pca_true[,2] <- silhouette(as.integer(metadata$CellType), dist = liger_dist_pca)[,3]

# load the workspace of MNN
PCs_MNN <- prcomp(corrected.mat, rank=pc_rank)$x
MNN_dist_pca <- dist(PCs_MNN)
sil_pca_true[,3] <- silhouette(as.integer(metadata$CellType), dist = MNN_dist_pca)[,3]

# load the workspace of scanorama
PCs_scanorama <- prcomp(scanorama_embedding, rank=pc_rank)$x
scanorama_dist_pca <- dist(PCs_scanorama)
sil_pca_true[,4] <- silhouette(as.integer(metadata$CellType), dist = scanorama_dist_pca)[,3]

# load the workspace of scVI
PCs_scVI <- prcomp(scVI_corrected, rank=pc_rank)$x
scVI_dist_pca <- dist(PCs_scVI)
sil_pca_true[,5] <- silhouette(as.integer(metadata$CellType), dist = scVI_dist_pca)[,3]

# load the workspace of seurat
PCs_Seurat <- Seurat_PCA[,1:10]
seurat_dist_pca <- dist(PCs_Seurat)
sil_pca_true[,6] <- silhouette(as.integer(metadata$CellType), dist = seurat_dist_pca)[,3]

# load the workspace of ZINBWaVE
PCs_ZINBWaVE <- prcomp(zinb_batch@W, rank=pc_rank)$x
ZINBWaVE_dist_pca <- dist(PCs_ZINBWaVE)
sil_pca_true[,7] <- silhouette(as.integer(metadata$CellType), dist = ZINBWaVE_dist_pca)[,3]

# Drawing the boxplot plot
colnames(sil_pca_true) <- c("BUSseq", "LIGER", "MNN", "Seurat", "Scanorama", "scVI", "ZINB-WaVE")
sil_com_pca_true<-melt(sil_pca_true)
colnames(sil_com_pca_true) <- c("No_Cell","Method","Silhouette_Coefficient")

jpeg("Image/boxplot_of_pancreas_Silhouette_coefs_10PCs.jpg",width = 900, height = 600,quality = 100)
p <- ggplot(sil_com_pca_true, aes(x=Method, y=Silhouette_Coefficient,fill=Method)) +
  geom_boxplot() + xlab(NULL) + ylab(NULL)+theme_bw() +
  theme(text = element_text(size=48),
        axis.text.x =element_text(face = "bold", #color = "#993333", 
                                  size = 36,angle = 90), 
        axis.text.y = element_text(face = "bold", #color = "#993333", 
                                   size = 36),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

write.table(sil_tSNE_true,"Results/Silhouette_coefs_10PCs.txt",row.names = F)
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

# Uncorrected #

set.seed(123)
all.dists.unc <- as.matrix(dist(log1p(t(data_raw))))
tsne_uncorrected <- Rtsne(all.dists.unc, is_distance=TRUE, perplexity = 30)

unc_by_celltype<- "Image/tsne_hemat_uncorrected_by_celltype.jpeg"
dat_frame <- data.frame(Var1 = tsne_uncorrected$Y[,1], 
                        Var2 = tsne_uncorrected$Y[,2], 
                        col = celltype_factor)

jpeg(unc_by_celltype,width = 800, height = 600, quality = 100)
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

unc_by_batch <- "Image/tsne_hemat_uncorrected_by_batch.jpeg"
dat_frame <- data.frame(Var1 = tsne_uncorrected$Y[,1], 
                        Var2 = tsne_uncorrected$Y[,2], 
                        col = batch_factor)

jpeg(unc_by_batch,width = 800, height = 600, quality = 100)
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

# Uncorrected #
pca.unc <- prcomp(log1p(t(data_raw)), rank=2)

unc_PCA<- "Image/pca_hemat_uncorrected.jpeg"
dat_frame <- data.frame(Var1 = pca.unc$x[,1], 
                        Var2 = pca.unc$x[,2], 
                        col = celltype_factor)

jpeg(unc_PCA,width = 800, height = 600, quality = 100)
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

# legend
jpeg(file="Image/legend_hemat_v1.jpeg", width=800, height=600,quality = 100)
plot(c(-5,5),c(-4,4),type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-5, y=4, legend=c("GSE72857","GSE81682"), pch=20, cex=2.5, col=color_by_batch, pt.cex = 6, title = "Batch", bty="n")
legend(x=-1, y=4, legend=c("LTHSC","MPP","LMPP","CMP","MEP","GMP","other"), pch=20, cex=2.5, pt.cex = 6, col=color_by_celltype, title = "Cell Type", bty="n")
dev.off()


# Store the workspace
save.image("comparison_workspace.RData")