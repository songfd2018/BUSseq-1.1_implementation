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
setwd("BatchEffectCorrectionComparison/liger/")

######################
# Load pancreas Data #
######################
# Loading pancreas count data
countdata <- read.table("../../RawCountData/count_data_pancreas_v1.txt") 

# Load dimension
dim <- unlist(read.table("../../RawCountData/dim_pancreas_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_pancreas_v1.txt")

# Load gene_list
gene_list <- unlist(read.table("../../RawCountData/gene_list_pancreas_v1.txt",stringsAsFactors = FALSE))

rownames(countdata) <- gene_list
colnames(countdata) <- rownames(metadata)
pancreas.list <- list(GSE81076 = countdata[,1:nb[1]], 
                      GSE85241 = countdata[,1:nb[2]+nb[1]], 
                      GSE86473 = countdata[,1:nb[3]+sum(nb[1:2])], 
                      EMTAB5061 = countdata[,1:nb[4]+sum(nb[1:3])])

################################################################################################################
# Conduct analysis referring to https://github.com/MacoskoLab/liger/blob/master/vignettes/walkthrough_pbmc.pdf #
################################################################################################################
# Create liger object
liger_pancreas <- createLiger(pancreas.list, remove.missing = FALSE)

# Normalization
liger_pancreas <- normalize(liger_pancreas)
liger_pancreas@var.genes <- gene_list
liger_pancreas <- scaleNotCenter(liger_pancreas)

# running suggestK on multiple cores can greatly decrease the runtime
# looking for the section of the plot where this metric stops increasing as sharply
k.suggest <- suggestK(liger_pancreas, num.cores = 5, gen.new = T, return.data = T, plot.log2 = F)


# Take the lowest objective of three factorizations with different initializations
# Multiple restarts are recommended for initial analyses since iNMF is non-deterministic
liger_pancreas <- optimizeALS(liger_pancreas, k=30, thresh = 5e-5, nrep = 3)

liger_pancreas <- runTSNE(liger_pancreas, use.raw = T)
p1 <- plotByDatasetAndCluster(liger_pancreas, return.plots = T)
# Plot by dataset
print(p1[[1]])

# To better integrate the datasets, we perform a quantile alignment step. This process first identifies similarly
# loading cells across datasets by building a similarity graph based on shared factor neighborhoods.
liger_pancreas <- quantileAlignSNF(liger_pancreas, resolution = 0.4, small.clust.thresh = 20)
w_liger <- liger_pancreas@alignment.clusters

#######
# ARI #
#######
ARI_liger <- adjustedRandIndex(metadata$CellType,w_liger)

kclust_liger <- pam(liger_pancreas@H.norm, k = 7)
w_liger_kmeans <- kclust_liger$clustering
ARI_liger_kmeans <- adjustedRandIndex(w_liger_kmeans,metadata$CellType)

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
loading_liger <- liger_pancreas@H.norm
liger_dist <- dist(loading_liger)
all.dists.liger <- as.matrix(liger_dist)
tsne_liger_dist <- Rtsne(all.dists.liger, is_distance=TRUE, perplexity = 30)

liger_by_celltype<- "Image/tsne_pancreas_liger_by_celltype.jpg"
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

liger_by_batch <- "Image/tsne_pancreas_liger_by_batch.jpg"
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

# Store the workspace
save(ARI_liger, ARI_liger_kmeans, tsne_liger_dist, loading_liger, file = "liger_results.RData")