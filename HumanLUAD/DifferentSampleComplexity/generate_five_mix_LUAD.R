rm(list = ls())
# run BUSseq on the benchmarked dataset
setwd("DifferentSampleComplexity")

library(scran)
library(scater)
library(ggplot2)
library(Rtsne)
library(RColorBrewer)
library(xtable)
################
# load dataset #
################
load("../sc_mixology-master/data/9cellmix_qc.RData")

###################################################
# feasure selection to generate the observed data #
###################################################
Num_hvg <- 6000

proj <- "LUAD"
ver <- 3
interested_type <- c("9-0-0","0-9-0","0-0-9","5-0-4","4-0-5")
# following the workflow of MNN
# as http://bioconductor.org/packages/3.7/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html



## cell_mix1 dataset
## cell selection
celltype <- paste0(sce_SC1_qc$HCC827, "-", sce_SC1_qc$H1975, "-", sce_SC1_qc$H2228)
sce_SC1_filter <- sce_SC1_qc[,celltype %in% interested_type ]

## normlization
clusters <- quickCluster(sce_SC1_filter, min.mean=0.1)
sce_SC1_filter <- computeSumFactors(sce_SC1_filter, min.mean=0.1, clusters=clusters)
sce_SC1_filter <- normalize(sce_SC1_filter)

var.fit <- trendVar(sce_SC1_filter, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_SC1_filter, var.fit)

plot(var.out$mean, var.out$total, xlab="Mean log-expression",
     ylab="Variance of log-expression", pch=16)

## select HVGs
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:Num_hvg], ]
hvg_mix1 <- rownames(hvg.out)

## cell_mix2 dataset
## cell selection
celltype <- paste0(sce_SC2_qc$HCC827, "-", sce_SC2_qc$H1975, "-", sce_SC2_qc$H2228)
sce_SC2_filter <- sce_SC2_qc[,celltype %in% interested_type ]

## normlization
clusters <- quickCluster(sce_SC2_filter, min.mean=0.1)
sce_SC2_filter <- computeSumFactors(sce_SC2_filter, min.mean=0.1, clusters=clusters)
sce_SC2_filter <- normalize(sce_SC2_filter)

var.fit <- trendVar(sce_SC2_filter, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_SC2_filter, var.fit)

plot(var.out$mean, var.out$total, xlab="Mean log-expression",
     ylab="Variance of log-expression", pch=16)

## select HVGs
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:Num_hvg], ]
hvg_mix2 <- rownames(hvg.out)

## cell_mix3 dataset
## cell selection
celltype <- paste0(sce_SC3_qc$HCC827, "-", sce_SC3_qc$H1975, "-", sce_SC3_qc$H2228)
sce_SC3_filter <- sce_SC3_qc[,celltype %in% interested_type ]

## normlization
clusters <- quickCluster(sce_SC3_filter, min.mean=0.1)
sce_SC3_filter <- computeSumFactors(sce_SC3_filter, min.mean=0.1, clusters=clusters)
sce_SC3_filter <- normalize(sce_SC3_filter)

var.fit <- trendVar(sce_SC3_filter, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_SC3_filter, var.fit)

plot(var.out$mean, var.out$total, xlab="Mean log-expression",
     ylab="Variance of log-expression", pch=16)

## select HVGs
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:Num_hvg], ]
hvg_mix3 <- rownames(hvg.out)

## cell_mix4 dataset
## cell selection
celltype <- paste0(sce_SC4_qc$HCC827, "-", sce_SC4_qc$H1975, "-", sce_SC4_qc$H2228)
sce_SC4_filter <- sce_SC4_qc[,celltype %in% interested_type ]

## normlization
clusters <- quickCluster(sce_SC4_filter, min.mean=0.1)
sce_SC4_filter <- computeSumFactors(sce_SC4_filter, min.mean=0.1, clusters=clusters)
sce_SC4_filter <- normalize(sce_SC4_filter)

var.fit <- trendVar(sce_SC4_filter, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_SC4_filter, var.fit)

plot(var.out$mean, var.out$total, xlab="Mean log-expression",
     ylab="Variance of log-expression", pch=16)

## select HVGs
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:Num_hvg], ]
hvg_mix4 <- rownames(hvg.out)


# Combined HVG
common.hvg <- intersect(hvg_mix1, intersect(hvg_mix2, intersect(hvg_mix3, hvg_mix4)))


####################
# Combined dataset #
####################
# index of common HVG in mix1 dataset
match_index <- match(common.hvg, rownames(sce_SC1_filter))
CountMat_mix1 <- counts(sce_SC1_filter)[match_index,]

# index of common HVG in mix2 dataset
match_index <- match(common.hvg, rownames(sce_SC2_filter))
CountMat_mix2 <- counts(sce_SC2_filter)[match_index,]

# index of common HVG in mix3 dataset
match_index <- match(common.hvg, rownames(sce_SC3_filter))
CountMat_mix3 <- counts(sce_SC3_filter)[match_index,]

# index of common HVG in mix4 dataset
match_index <- match(common.hvg, rownames(sce_SC4_filter))
CountMat_mix4 <- counts(sce_SC4_filter)[match_index,]

nb <- c(ncol(CountMat_mix1), ncol(CountMat_mix2), ncol(CountMat_mix3), ncol(CountMat_mix4))

LUADmixCounts <- cbind(CountMat_mix1, CountMat_mix2, CountMat_mix3, CountMat_mix4)
colnames(LUADmixCounts) <- c(paste0("mix1_",colnames(CountMat_mix1)),
                             paste0("mix2_",colnames(CountMat_mix2)),
                             paste0("mix3_",colnames(CountMat_mix3)),
                             paste0("mix4_",colnames(CountMat_mix4)))

# metadata
Num_Cell1 <- c(sce_SC1_filter$HCC827, sce_SC2_filter$HCC827, sce_SC3_filter$HCC827, sce_SC4_filter$HCC827)
Num_Cell2 <- c(sce_SC1_filter$H1975, sce_SC2_filter$H1975, sce_SC3_filter$H1975, sce_SC4_filter$H1975)
Num_Cell3 <- c(sce_SC1_filter$H2228, sce_SC2_filter$H2228, sce_SC3_filter$H2228, sce_SC4_filter$H2228)

batch_label <- rep(c("Mix1","Mix2","Mix3","Mix4"), nb)

metadata <- data.frame(Batch = batch_label, 
                       Num_Cell1 = Num_Cell1,
                       Num_Cell2 = Num_Cell2,
                       Num_Cell3 = Num_Cell3)
rownames(metadata) <- colnames(LUADmixCounts)


#store count data matrix
if(!dir.exists("RawCountData")){
  dir.create("RawCountData")
}
write.table(LUADmixCounts,file=paste0("RawCountData/count_data_",proj,"_v",ver,".txt"),row.names = FALSE, col.names = FALSE)
write.table(metadata,file=paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

N <- ncol(LUADmixCounts)
G <- nrow(LUADmixCounts)
B <- length(nb)
write.table(c(N,G,B,nb),file=paste0("RawCountData/dim_",proj,"_v",ver,".txt"),row.names = FALSE, col.names = FALSE)
write.table(common.hvg,paste0("RawCountData/gene_list_",proj,"_v",ver,".txt"),row.names = FALSE, col.names = FALSE)

save.image(file = paste0("CountData_",proj,"_v",ver,".RData"))

if(!dir.exists("five_mix")){
  dir.create("five_mix")
}


# Draw tSNE and PCA plot for the raw count data
if(!dir.exists("Image")){
  dir.create("Image")
}

color_by_celltype <- c(brewer.pal(3,"Set3"),"black")
celltypes <- c("HCC827", "H1975", "H2228", "Equal")
celltype_factor <- rep(NA, N)
for(i in 1:N){
  if(metadata[i,2] == metadata[i,3] & metadata[i,2] == metadata[i,4] & metadata[i,3] == metadata[i,4]){
    celltype_factor[i] <- "Equal"
  }else{
    cell_index <- which.max(metadata[i,2:4])
    celltype_factor[i] <- celltypes[cell_index]
  }
}
celltype_factor <- factor(celltype_factor, levels = celltypes)

#####set batch colorings
color_by_batch <- c("#EB4334","#FBBD06","#35AA53","#4586F3")
batch_factor <- factor(rep(paste0("mix",1:B),nb))


# Uncorrected #
set.seed(123)
all.dists.unc <- as.matrix(dist(log1p(t(LUADmixCounts))))
tsne_uncorrected <- Rtsne(all.dists.unc, is_distance=TRUE, perplexity = 30)

uncorrected_by_celltype<- paste0("Image/tsne_",proj,"_uncorrected_by_celltype_v",ver,".jpg")
dat_frame <- data.frame(Var1 = tsne_uncorrected$Y[,1], 
                        Var2 = tsne_uncorrected$Y[,2], 
                        col = celltype_factor)

jpeg(uncorrected_by_celltype,width = 800, height = 600, quality = 100)
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
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"))
dev.off()

uncorrected_by_batch <- paste0("Image/tsne_",proj,"_uncorrected_by_batch_v",ver,".jpg")
dat_frame <- data.frame(Var1 = tsne_uncorrected$Y[,1], 
                        Var2 = tsne_uncorrected$Y[,2], 
                        col = batch_factor)

jpeg(uncorrected_by_batch,width = 800, height = 600, quality = 100)
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
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"))
dev.off()

##################
# Draw PCA plots #
##################
# Uncorrected #
pca.unc <- prcomp(log1p(t(LUADmixCounts)), rank=2)
unc_PCA_by_celltype<- paste0("Image/pca_",proj,"_uncorrected_by_celltype_v",ver,".jpg")
dat_frame <- data.frame(Var1 = pca.unc$x[,1], 
                        Var2 = pca.unc$x[,2], 
                        col = celltype_factor)
jpeg(unc_PCA_by_celltype,width = 800, height = 600, quality = 100)
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
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"))
dev.off()

unc_PCA_by_batch<- paste0("Image/pca_",proj,"_uncorrected_by_batch_v",ver,".jpg")
dat_frame <- data.frame(Var1 = pca.unc$x[,1], 
                        Var2 = pca.unc$x[,2], 
                        col = batch_factor)
jpeg(unc_PCA_by_batch,width = 800, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_batch,0.6)) +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"))
dev.off()
