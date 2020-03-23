rm(list=ls())
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)
library(reshape2)

# coloring
library(scales)
require(WGCNA)
library(RColorBrewer)

###################
# Load results #
###################
# Working directory
setwd("MainSimulation")

load("BUSseq/BUSseq_results.RData")
comparison_list <- c("liger", "MNN", "scanorama", "scVI", "Seurat", "ZINBWaVE")
num_comparison <- length(comparison_list)

for(m in 1:num_comparison){
  load(paste0("Comparison/",comparison_list[m],"/",comparison_list[m],"_results.RData"))
}

# load the dimentsion information
dim <- unlist(read.table("RawCountData/dim_simulation_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[1:B + 3]

# Load metadata
metadata <- read.table("RawCountData/metadata_simulation_v1.txt")

# Load gene_list
gene_list <- paste0("gene-",1:G)

# Load raw count data 
data_raw <- as.matrix(read.table("RawCountData/count_data_simulation_v1.txt"))

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

###############################################################################
# Silhouette coefficient of tSNE by true cell type labels of each method #
###############################################################################
# Storing the Silhouette coefficients
sil_tSNE_true <- matrix(NA,N,num_comparison+1)
colnames(sil_tSNE_true) <- names(ARI_values)

BUSseq_dist_tSNE <- dist(tsne_BUSseq_dist$Y)
sil_tSNE_true[,1] <- silhouette(as.integer(metadata$celltype), dist = BUSseq_dist_tSNE)[,3]

liger_dist_tSNE <- dist(tsne_liger_dist$Y)
sil_tSNE_true[,2] <- silhouette(as.integer(metadata$celltype), dist = liger_dist_tSNE)[,3]

MNN_dist_tSNE <- dist(tsne_MNN_dist$Y)
sil_tSNE_true[,3] <- silhouette(as.integer(metadata$celltype), dist = MNN_dist_tSNE)[,3]

scanorama_dist_tSNE <- dist(tsne_scanorama_dist$Y)
sil_tSNE_true[,4] <- silhouette(as.integer(metadata$celltype), dist = scanorama_dist_tSNE)[,3]

scVI_dist_tSNE <- dist(tsne_scVI_dist$Y)
sil_tSNE_true[,5] <- silhouette(as.integer(metadata$celltype), dist = scVI_dist_tSNE)[,3]

Seurat_dist_tSNE <- dist(tsne_Seurat_dist$Y)
sil_tSNE_true[,6] <- silhouette(as.integer(metadata$celltype), dist = Seurat_dist_tSNE)[,3]

ZINBWaVE_dist_tSNE <- dist(tsne_ZINBW_dist$Y)
sil_tSNE_true[,7] <- silhouette(as.integer(metadata$celltype), dist = ZINBWaVE_dist_tSNE)[,3]

write.table(sil_tSNE_true,"Results/Silhouette_coefs.txt",row.names = F)

# Drawing the boxplot plot
colnames(sil_tSNE_true) <- c("BUSseq", "LIGER", "MNN", "Seurat", "Scanorama", "scVI", "ZINB-WaVE")
sil_com_tSNE_true<-melt(sil_tSNE_true)
colnames(sil_com_tSNE_true) <- c("No_Cell","Method","Silhouette_Coefficient")

jpeg("Image/boxplot_of_Silhouette_coefs.jpg",width = 900, height = 600,quality = 100)
p <- ggplot(sil_com_tSNE_true, aes(x=Method, y=Silhouette_Coefficient,fill=Method)) +
  geom_boxplot() + xlab(NULL) + ylab(NULL)+theme_bw() +
  theme(text = element_text(size=48),
        axis.text.x =element_text(face = "bold",
                                  size = 36,angle = 90), 
        axis.text.y = element_text(face = "bold",
                                   size = 36),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

#################################
# scatter plots (t-SNE and PCA) #
#################################
#set cell type colorings
# batch color bar
color_by_batch<-c("#EB4334","#FBBD06","#35AA53","#4586F3")
batch_factor <- factor(paste0("Batch ",rep(1:B,nb)))

# cell type color bar
color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978","#8097D3","#9C7ACE")
celltype_factor <- factor(paste0("Celltype ",metadata$celltype))

# Uncorrected #
set.seed(123)
all.dists.unc <- as.matrix(dist(log1p(t(data_raw))))
tsne_uncorrected <- Rtsne(all.dists.unc, is_distance=TRUE, perplexity = 30)

unc_by_celltype<- "Image/tsne_simulation_uncorrected_by_celltype.jpeg"
dat_frame <- data.frame(Var1 = tsne_uncorrected$Y[,1], 
                        Var2 = tsne_uncorrected$Y[,2], 
                        col = celltype_factor)

jpeg(unc_by_celltype,width = 900, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_celltype,0.6)) +
  xlab("tSNE 1") + ylab("tSNE 2") +
  theme(axis.text.x = element_text(face = "bold",
                                   size = 44),
        axis.text.y = element_text(face = "bold",
                                   size = 44),
        axis.title=element_text(size=48,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.position = "none")
dev.off()

unc_by_batch <- "Image/tsne_simulation_uncorrected_by_batch.jpeg"
dat_frame <- data.frame(Var1 = tsne_uncorrected$Y[,1], 
                        Var2 = tsne_uncorrected$Y[,2], 
                        col = batch_factor)

jpeg(unc_by_batch,width = 900, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_batch,0.6)) +
  xlab("tSNE 1") + ylab("tSNE 2") +
  theme(axis.text.x = element_text(face = "bold",
                                   size = 36),
        axis.text.y = element_text(face = "bold", 
                                   size = 36),
        axis.title=element_text(size=44,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.position = "none")
dev.off()

##################
# Draw PCA plots #
##################

# Uncorrected #
pca.unc <- prcomp(log1p(t(data_raw)), rank=2)
unc_PCA<- "Image/pca_simulation_uncorrected.jpeg"
dat_frame <- data.frame(Var1 = pca.unc$x[,1], 
                        Var2 = pca.unc$x[,2], 
                        col = celltype_factor)
jpeg(unc_PCA,width = 900, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_celltype,0.6)) +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold",
                                   size = 36),
        axis.text.y = element_text(face = "bold", 
                                   size = 36),
        axis.title=element_text(size=44,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.position = "none")
dev.off()

# legend
K <- 5
jpeg("Image/legend_simulation_v1.jpeg",width = 600, height = 400, quality = 100)
plot(c(-5,5),c(-4,4),type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-5, y=4, legend=paste0("Batch ",1:B), pch=20, cex=2.5, col=color_by_batch, title = "Batch", bty="n")
legend(x=-0, y=4, legend=paste0("Celltype ",1:K), pch=20, cex=2.5, col=color_by_celltype, title = "Cell Type", bty="n")
dev.off()

# Store the workspace
save.image("comparison_workspace.RData")