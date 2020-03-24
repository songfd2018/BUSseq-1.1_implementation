# library(devtools)
# install_github('MacoskoLab/liger')
rm(list=ls())
library(liger)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)

# Working directory
setwd("BatchEffectCorrectionComparison/liger/")

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

# Load gene_list
gene_list <- unlist(read.table("../../RawCountData/gene_list_hemat_v1.txt",stringsAsFactors = FALSE))

rownames(countdata) <- gene_list
colnames(countdata) <- paste0("Batch",rep(1:2,nb),"_Cell",c(1:nb[1],1:nb[2]))
hemat.list <- list(GSE72857 = countdata[,1:nb[1]], GSE81682 = countdata[,1:nb[2]+nb[1]])


################################################################################################################
# Conduct analysis referring to https://github.com/MacoskoLab/liger/blob/master/vignettes/walkthrough_pbmc.pdf #
################################################################################################################
# Create liger object
liger_hemat <- createLiger(hemat.list, remove.missing = FALSE)

# Normalization
liger_hemat <- normalize(liger_hemat)
liger_hemat@var.genes <- gene_list
liger_hemat <- scaleNotCenter(liger_hemat)

# running suggestK on multiple cores can greatly decrease the runtime
# looking for the section of the plot where this metric stops increasing as sharply
# k.suggest <- suggestK(liger_hemat, num.cores = 5, gen.new = T, return.data = T, plot.log2 = F)

# Take the lowest objective of three factorizations with different initializations
# Multiple restarts are recommended for initial analyses since iNMF is non-deterministic
liger_hemat <- optimizeALS(liger_hemat, k=22, thresh = 5e-5, nrep = 3)

liger_hemat <- runTSNE(liger_hemat, use.raw = T)
p1 <- plotByDatasetAndCluster(liger_hemat, return.plots = T)
# Plot by dataset
print(p1[[1]])

# To better integrate the datasets, we perform a quantile alignment step. This process first identifies similarly
# loading cells across datasets by building a similarity graph based on shared factor neighborhoods.
liger_hemat <- quantileAlignSNF(liger_hemat, resolution = 0.4, small.clust.thresh = 20)

metadata <- read.table(file = "../../RawCountData/hemat_metadata.txt")
GSE72857_celltype <- droplevels(metadata$CellType[1:nb[1]])
GSE81682_celltype <- droplevels(metadata$CellType[1:nb[2] + nb[1]])


print(p_a[[2]])

# This function uses the V loadings to identify dataset-specific genes and a
# combination of the W and V loadings to identify shared markers
# identify some shared and dataset-specific markers for each factor and plot them to help in cluster
# annotation.
markers <- getFactorMarkers(liger_hemat, dataset1='GSE72857', dataset2='GSE81682',
                            num.genes = 10)
marker_genes <- markers$shared[order(markers$shared$p_value),]
w_liger <- liger_hemat@alignment.clusters

#######
# ARI #
#######
ARI_liger <- adjustedRandIndex(metadata$CellType,w_liger)

kclust_liger <- pam(liger_hemat@H.norm, k = 7)
w_liger_kmeans <- kclust_liger$clustering
ARI_liger_kmeans <- adjustedRandIndex(w_liger_kmeans,metadata$CellType)

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
loading_liger <- liger_hemat@H.norm
liger_dist <- dist(loading_liger)
all.dists.liger <- as.matrix(liger_dist)
tsne_liger_dist <- Rtsne(all.dists.liger, is_distance=TRUE, perplexity = 30)

liger_by_celltype<- "Image/tsne_hemat_liger_by_celltype.jpg"
dat_frame <- data.frame(Var1 = tsne_liger_dist$Y[,1], 
                        Var2 = tsne_liger_dist$Y[,2], 
                        col = celltype_factor)

jpeg(liger_by_celltype,width = 800, height = 600, quality = 100)
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

liger_by_batch <- "Image/tsne_hemat_liger_by_batch.jpg"
dat_frame <- data.frame(Var1 = tsne_liger_dist$Y[,1], 
                        Var2 = tsne_liger_dist$Y[,2], 
                        col = batch_factor)

jpeg(liger_by_batch,width = 800, height = 600, quality = 100)
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

liger_PCA<- "Image/pca_hemat_liger.jpeg"
dat_frame <- data.frame(Var1 = pca.liger$x[,1], 
                        Var2 = pca.liger$x[,2], 
                        col = celltype_factor)

jpeg(liger_PCA,width = 800, height = 600, quality = 100)
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
save(ARI_liger, ARI_liger_kmeans, tsne_liger_dist, loading_liger, file = "liger_results.RData")