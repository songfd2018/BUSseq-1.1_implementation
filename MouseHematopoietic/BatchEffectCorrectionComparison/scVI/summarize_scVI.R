# Apply scVI to the hematopoietic study.
rm(list=ls())
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)

set.seed(123)

# Working directory
setwd("BatchEffectCorrectionComparison/scVI/")

###################
# Load hemat Data #
###################
# load the dimentsion information
dim <- unlist(read.table("../../RawCountData/dim_hemat_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[1:B + 3]

# Loading the true cell type indicators
metadata <- read.table("../../RawCountData/metadata_hemat_v1.txt")

##############################
# load the inference by scVI #
##############################
w_scVI <- unlist(read.table("scVI_hemat_v1_clusters_0716.txt"))
scVI_corrected <- read.table("scVI_hemat_v1_latent_0716.txt")

#######
# ARI #
#######
ARI_scVI <- adjustedRandIndex(metadata$CellType, w_scVI)

kclust_scVI <- pam(scVI_corrected, k = 7)
w_scVI_kmeans <- kclust_scVI$clustering
ARI_scVI_kmeans <- adjustedRandIndex(w_scVI_kmeans,metadata$CellType)

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
scVI_dist <- dist(scVI_corrected) 
all.dists.scVI <- as.matrix(scVI_dist)
tsne_scVI_dist <- Rtsne(all.dists.scVI, is_distance=TRUE, perplexity = 30)

scVI_by_celltype<- "Image/tsne_hemat_scVI_by_celltype.jpg"
dat_frame <- data.frame(Var1 = tsne_scVI_dist$Y[,1], 
                        Var2 = tsne_scVI_dist$Y[,2], 
                        col = celltype_factor)

jpeg(scVI_by_celltype,width = 800, height = 600, quality = 100)
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

scVI_by_batch <- "Image/tsne_hemat_scVI_by_batch.jpg"
dat_frame <- data.frame(Var1 = tsne_scVI_dist$Y[,1], 
                        Var2 = tsne_scVI_dist$Y[,2], 
                        col = batch_factor)

jpeg(scVI_by_batch,width = 800, height = 600, quality = 100)
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

scVI_PCA<- "Image/pca_hemat_scVI.jpeg"
dat_frame <- data.frame(Var1 = pca.scVI$x[,1], 
                        Var2 = pca.scVI$x[,2], 
                        col = celltype_factor)

jpeg(scVI_PCA,width = 800, height = 600, quality = 100)
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
save(ARI_scVI, ARI_scVI_kmeans, tsne_scVI_dist, scVI_corrected, file = "scVI_results.RData")