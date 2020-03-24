#Apply MNN to the hematopoietic study.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("scran")

rm(list=ls())
library(scran)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)

set.seed(123)
# Working directory
setwd("BatchEffectCorrectionComparison/MNN/")

###################
# Load hemat Data #
###################
# Loading the file name list of all hemat count data
countdata <- read.table("../../RawCountData/count_data_hemat_v1.txt") 

dim <- unlist(read.table("../../RawCountData/dim_hemat_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_hemat_v1.txt")
gene_list <- unlist(read.table("../../RawCountData/gene_list_hemat_v1.txt",stringsAsFactors = FALSE))

colnames(countdata) <- rownames(metadata)
rownames(countdata) <- gene_list

HematCounts <- list()
index <- 0
for(b in 1:B){
  HematCounts[[b]] <- as.matrix(countdata[,index + 1:nb[b]])
  index <- index + nb[b]
}

#######################################
# Apply MNN to the Hematopoietic Data #
#######################################
data_MNN <- HematCounts
# Referring to https://github.com/MarioniLab/MNN2017/blob/master/Haematopoiesis/plotCorrections.R
# Normalization
data_MNN_normalized<-list()
for(b in 1:B){
  high.abF <- scater::calcAverage(data_MNN[[b]]) > 1
  clustF <- quickCluster(data_MNN[[b]], min.size=10 , method="igraph", subset.row=high.abF)
  sizeF <- computeSumFactors(data_MNN[[b]], sizes=seq(11, 81, 5), cluster=clustF, subset.row=high.abF)
  data_MNN_normalized[[b]] <- t(t(data_MNN[[b]])/sizeF)
}

# Rescaling the first dataset to match the coverage of the second.
ave<-list()
for(b in 1:B){
  ave[[b]] <- rowMeans(data_MNN_normalized[[b]])
  if(b>1){
    data_MNN_normalized[[b]] <- data_MNN_normalized[[b]] * median(ave[[1]]/ave[[b]])
  }
}

# Performing log-transformation and save results to file.
log_data_MNN<-list()
for(b in 1:B){
  log_data_MNN[[b]] <- log1p(data_MNN_normalized[[b]])
}

# Performing the correction with MNN (turned down the sigma to improve mixing).

mnn.out<-mnnCorrect(log_data_MNN[[1]], log_data_MNN[[2]],k=20, sigma=0.1,cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE,compute.angle=TRUE)

X.mnn <- cbind(mnn.out$corrected[[1]], mnn.out$corrected[[2]])
t.mnn <- t(X.mnn)

##############
# Clustering #
##############
# Referring to http://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/batch.html
omat <- do.call(cbind, log_data_MNN)
sce <- SingleCellExperiment(list(logcounts=omat))
reducedDim(sce, "MNN") <- t.mnn
sce$Batch <- rep(paste0("Batch",1:2),nb)

start_time <- Sys.time()
snn.gr <- buildSNNGraph(sce, use.dimred="MNN")
clusters <- igraph::cluster_walktrap(snn.gr)
end_time <- Sys.time()
time_consumption <- end_time - start_time

table(clusters$membership, sce$Batch)

# The estimated cell type indicators by MNN
w_MNN <- factor(clusters$membership)

#######
# ARI #
#######
ARI_MNN <- adjustedRandIndex(metadata$CellType,w_MNN)

kclust_MNN <- pam(t(omat), k = 7)
w_MNN_kmeans <- kclust_MNN$clustering
ARI_MNN_kmeans <- adjustedRandIndex(w_MNN_kmeans,metadata$CellType)

#######################################################################
# Silhouette coefficient by estimated cell type labels of each method #
#######################################################################
MNN_dist <- dist(t.mnn)
sil_MNN <- silhouette(as.integer(w_MNN), dist = MNN_dist)
sil_MNN_true <- silhouette(as.integer(metadata$CellType), dist = MNN_dist)

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
all.dists.mnn <- as.matrix(MNN_dist)
tsne_MNN_dist <- Rtsne(all.dists.mnn, is_distance=TRUE, perplexity = 30)

MNN_by_celltype<- "Image/tsne_hemat_MNN_by_celltype.jpg"
dat_frame <- data.frame(Var1 = tsne_MNN_dist$Y[,1], 
                        Var2 = tsne_MNN_dist$Y[,2], 
                        col = celltype_factor)

jpeg(MNN_by_celltype,width = 800, height = 600, quality = 100)
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

MNN_by_batch <- "Image/tsne_hemat_MNN_by_batch.jpg"
dat_frame <- data.frame(Var1 = tsne_MNN_dist$Y[,1], 
                        Var2 = tsne_MNN_dist$Y[,2], 
                        col = batch_factor)

jpeg(MNN_by_batch,width = 800, height = 600, quality = 100)
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
pca.mnn <- prcomp(t.mnn, rank=2)

MNN_PCA<- "Image/pca_hemat_MNN.jpeg"
dat_frame <- data.frame(Var1 = pca.mnn$x[,1], 
                        Var2 = pca.mnn$x[,2], 
                        col = celltype_factor)

jpeg(MNN_PCA,width = 800, height = 600, quality = 100)
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
save(ARI_MNN, ARI_MNN_kmeans, tsne_MNN_dist, t.mnn, file = "MNN_results.RData")