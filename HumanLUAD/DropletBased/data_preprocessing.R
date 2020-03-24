rm(list = ls())
# run BUSseq on the benchmarked dataset
setwd("DropletBased")

library(scran)
library(scater)
library(RColorBrewer)
################
# load dataset #
################
load("../sc_mixology-master/data/sincell_with_class.RData")

###################################################
# feasure selection to generate the observed data #
###################################################
Num_hvg <- 6000

proj <- "LUAD"
ver <- 1

# following the workflow of MNN
# as http://bioconductor.org/packages/3.7/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html


## 10x dataset
clusters <- quickCluster(sce_sc_10x_qc, min.mean=0.1)
sce_sc_10x_qc <- computeSumFactors(sce_sc_10x_qc, min.mean=0.1, clusters=clusters)
sce_sc_10x_qc <- normalize(sce_sc_10x_qc)

var.fit <- trendVar(sce_sc_10x_qc, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_sc_10x_qc, var.fit)

plot(var.out$mean, var.out$total, xlab="Mean log-expression", 
     ylab="Variance of log-expression", pch=16)

hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:Num_hvg], ]
hvg_10x <- rownames(hvg.out)

count_hvg <- counts(sce_sc_10x_qc)[rownames(sce_sc_10x_qc) %in% hvg_10x,]
pca_10x <- prcomp(log1p(t(count_hvg)), rank=2)
file_10x_PCA<- paste0("Image/pca_",proj,"_10x_v",ver,".jpeg")
dat_frame <- data.frame(Var1 = pca_10x$x[,1], 
                        Var2 = pca_10x$x[,2], 
                        col = factor(sce_sc_10x_qc$cell_line_demuxlet))

jpeg(file_10x_PCA,width = 800, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  # scale_colour_manual(values = alpha(color_by_celltype,0.6)) +
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

## CELseq dataset
clusters <- quickCluster(sce_sc_CELseq2_qc, min.mean=0.1)
sce_sc_CELseq2_qc <- computeSumFactors(sce_sc_CELseq2_qc, min.mean=0.1, clusters=clusters)
sce_sc_CELseq2_qc <- normalize(sce_sc_CELseq2_qc)

var.fit <- trendVar(sce_sc_CELseq2_qc, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_sc_CELseq2_qc, var.fit)

plot(var.out$mean, var.out$total, xlab="Mean log-expression", 
     ylab="Variance of log-expression", pch=16)

hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:Num_hvg], ]
hvg_CELseq2 <- rownames(hvg.out)

count_hvg <- counts(sce_sc_CELseq2_qc)[rownames(sce_sc_CELseq2_qc) %in% hvg_CELseq2,]
pca_CELseq2 <- prcomp(log1p(t(count_hvg)), rank=2)
file_CELseq2_PCA<- paste0("Image/pca_",proj,"_CELseq2_v",ver,".jpeg")
dat_frame <- data.frame(Var1 = pca_CELseq2$x[,1], 
                        Var2 = pca_CELseq2$x[,2], 
                        col = factor(sce_sc_CELseq2_qc$cell_line_demuxlet))

jpeg(file_CELseq2_PCA,width = 800, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  # scale_colour_manual(values = alpha(color_by_celltype,0.6)) +
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

## Dropseq dataset
clusters <- quickCluster(sce_sc_Dropseq_qc, min.mean=0.1)
sce_sc_Dropseq_qc <- computeSumFactors(sce_sc_Dropseq_qc, min.mean=0.1, clusters=clusters)
sce_sc_Dropseq_qc <- normalize(sce_sc_Dropseq_qc)

var.fit <- trendVar(sce_sc_Dropseq_qc, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_sc_Dropseq_qc, var.fit)

plot(var.out$mean, var.out$total, xlab="Mean log-expression", 
     ylab="Variance of log-expression", pch=16)

hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:Num_hvg], ]
hvg_Dropseq <- rownames(hvg.out)

# Combined HVG
common.hvg <- intersect(hvg_10x, intersect(hvg_CELseq2, hvg_Dropseq))

####################
# Combined dataset #
####################
# index of common HVG in 10x dataset
match_index <- match(common.hvg, rownames(sce_sc_10x_qc))
CountMat_10x <- counts(sce_sc_10x_qc)[match_index,]

# index of common HVG in CELseq2 dataset
match_index <- match(common.hvg, rownames(sce_sc_CELseq2_qc))
CountMat_CELseq2 <- counts(sce_sc_CELseq2_qc)[match_index,]

# index of common HVG in Dropseq dataset
match_index <- match(common.hvg, rownames(sce_sc_Dropseq_qc))
CountMat_Dropseq <- counts(sce_sc_Dropseq_qc)[match_index,]


nb <- c(ncol(CountMat_CELseq2), ncol(CountMat_10x), ncol(CountMat_Dropseq))

LUADCounts <- cbind(CountMat_CELseq2, CountMat_10x, CountMat_Dropseq)
colnames(LUADCounts) <- c(paste0("CELseq2_",colnames(CountMat_CELseq2)),
                          paste0("10x_",colnames(CountMat_10x)),
                          paste0("Dropseq_",colnames(CountMat_Dropseq)))

# metadata
demuxlet_label <- c(sce_sc_CELseq2_qc$cell_line_demuxlet, sce_sc_10x_qc$cell_line_demuxlet, sce_sc_Dropseq_qc$cell_line_demuxlet)
prot_label <- rep(c("CELseq2","10x","Dropseq"), nb)

metadata <- data.frame(Protocol = prot_label, 
                       CellType = demuxlet_label)
rownames(metadata) <- colnames(LUADCounts)


#store count data matrix
if(!dir.exists("RawCountData")){
  dir.create("RawCountData")
}
write.table(LUADCounts,file=paste0("RawCountData/count_data_",proj,"_v",ver,".txt"),row.names = FALSE, col.names = FALSE)
write.table(metadata,file=paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

N <- ncol(LUADCounts)
G <- nrow(LUADCounts)
B <- length(nb)
write.table(c(N,G,B,nb),file=paste0("RawCountData/dim_",proj,"_v",ver,".txt"),row.names = FALSE, col.names = FALSE)
write.table(common.hvg,paste0("RawCountData/gene_list_",proj,"_v",ver,".txt"),row.names = FALSE, col.names = FALSE)

save.image(file = paste0("CountData_",proj,"_v",ver,".RData"))

# Draw tSNE and PCA plot for the raw count data
if(!dir.exists("Image")){
  dir.create("Image")
}

#####set cell type colorings
color_by_celltype <- c(brewer.pal(3,"Set3"))
celltype_factor <- metadata$CellType

#####set batch colorings
color_by_batch <- c("#EB4334","#FBBD06","#35AA53")
batch_factor <- factor(metadata$Protocol)

# Uncorrected #
set.seed(123)
all.dists.unc <- as.matrix(dist(log1p(t(LUADCounts))))
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
        legend.position = "none")
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
        legend.position = "none")
dev.off()

##################
# Draw PCA plots #
##################
# Uncorrected #
pca.unc <- prcomp(log1p(t(LUADCounts)), rank=2)
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
        legend.position = "none")
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
        legend.position = "none")
dev.off()

# legend
jpeg(paste0("Image/legend_",proj,"_v",ver,".jpg"),width = 800, height = 600, quality = 100)
plot(c(-8,8),c(-6,6),type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-8, y=6, legend=c("HCC827", "H1975", "H2228"), pch=20, cex=2.5, pt.cex = 6, col=color_by_celltype , title = "Cell Type", bty="n")
legend(x=0, y=6, legend=c("CELseq2","10x","Dropseq"), pch=20, cex=2.5, pt.cex = 6, col = color_by_batch, title = "Batch", bty="n")
dev.off()