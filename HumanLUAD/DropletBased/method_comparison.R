##############################################
# Apply benchmarked methods to the LUAD data #
##############################################
setwd("DropletBased")

################################
# Apply LIGER to the LUAD Data #
################################
install.packages('devtools')
library(devtools)
install_github('MacoskoLab/liger')

rm(list=ls())
library(liger)
library(mclust) # For ARI
library(cluster) # For Silhouette
set.seed(12345)
method <- "LIGER"
proj <- "LUAD"
ver <- 1

##################
# Load LUAD Data #
##################
# Loading LUAD count data
countdata <- read.table(paste0("RawCountData/count_data_",proj,"_v",ver,".txt")) 

# Load dimension
dim <- unlist(read.table(paste0("RawCountData/dim_",proj,"_v",ver,".txt")))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

# Load gene_list
gene_list <- unlist(read.table(paste0("RawCountData/gene_list_",proj,"_v",ver,".txt"),stringsAsFactors = FALSE))

rownames(countdata) <- gene_list
colnames(countdata) <- rownames(metadata)
LUAD.list <- list(TenX = countdata[,1:nb[1]], CELseq2 = countdata[,1:nb[2]+nb[1]], Dropseq = countdata[,1:nb[3]+sum(nb[1:2])])

################################################################################################################
# Conduct analysis referring to https://github.com/MacoskoLab/liger/blob/master/vignettes/walkthrough_pbmc.pdf #
################################################################################################################
# Create liger object
liger_LUAD <- createLiger(LUAD.list, remove.missing = FALSE)

# Normalization
liger_LUAD <- normalize(liger_LUAD)
liger_LUAD@var.genes <- gene_list
liger_LUAD <- scaleNotCenter(liger_LUAD)

# running suggestK on multiple cores can greatly decrease the runtime
# looking for the section of the plot where this metric stops increasing as sharply
k.suggest <- suggestK(liger_LUAD, num.cores = 5, gen.new = T, return.data = T, plot.log2 = F)


# Take the lowest objective of three factorizations with different initializations
# Multiple restarts are recommended for initial analyses since iNMF is non-deterministic
liger_LUAD <- optimizeALS(liger_LUAD, k=40, thresh = 5e-5, nrep = 3)

liger_LUAD <- runTSNE(liger_LUAD, use.raw = T)
p1 <- plotByDatasetAndCluster(liger_LUAD, return.plots = T)
# Plot by dataset
print(p1[[1]])

# To better integrate the datasets, we perform a quantile alignment step. This process first identifies similarly
# loading cells across datasets by building a similarity graph based on shared factor neighborhoods.
liger_LUAD <- quantileAlignSNF(liger_LUAD, resolution = 0.4, small.clust.thresh = 20)
w_liger <- liger_LUAD@alignment.clusters

#######
# ARI #
#######
ARI_liger <- adjustedRandIndex(metadata$CellType,w_liger)

save.image(paste0("Comparison/",method, "_workspace.RData"))

##############################
# Apply MNN to the LUAD Data #
##############################
rm(list=ls())
library(scran)
library(batchelor)
library(mclust) # For ARI
library(cluster) # For Silhouette
set.seed(12345)
method <- "MNN"
proj <- "LUAD"
ver <- 1

##################
# Load LUAD Data #
##################
# Loading LUAD count data
countdata <- read.table(paste0("RawCountData/count_data_",proj,"_v",ver,".txt")) 

# Load dimension
dim <- unlist(read.table(paste0("RawCountData/dim_",proj,"_v",ver,".txt")))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

# Load gene_list
gene_list <- unlist(read.table(paste0("RawCountData/gene_list_",proj,"_v",ver,".txt"),stringsAsFactors = FALSE))

rownames(countdata) <- gene_list
colnames(countdata) <- rownames(metadata)

data_MNN <- list()
index <- 0
for(b in 1:B){
  data_MNN[[b]] <- as.matrix(countdata[,index + 1:nb[b]])
  index <- index + nb[b]
}

# Referring to https://github.com/MarioniLab/MNN2017/blob/master/Pancreas
# Referring to http://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/batch.html
# Normalization
data_MNN_normalized<-list()

for(b in 1:B){
  high.abF <- scater::calcAverage(data_MNN[[b]]) > 1
  clustF <- quickCluster(data_MNN[[b]], min.size=10 , method="igraph", subset.row=high.abF)
  sizeF <- computeSumFactors(data_MNN[[b]], sizes=seq(11, 81, 5), cluster=clustF, subset.row=high.abF)
  data_MNN_normalized[[b]] <- t(t(data_MNN[[b]])/sizeF)
}

# MNN batch correction
Xmnn <- mnnCorrect(data_MNN_normalized[[1]],
                   data_MNN_normalized[[2]],
                   data_MNN_normalized[[3]],
                   svd.dim=0,
                   cos.norm.in=TRUE, cos.norm.out=TRUE,
                   var.adj=TRUE, 
                   k=20, sigma=0.1)

##############
# Clustering #
##############

omat <- do.call(cbind, data_MNN_normalized)
sce <- SingleCellExperiment(list(logcounts=log1p(omat)))
reducedDim(sce, "Corrected") <- t(assay(Xmnn))
sce$Batch <- rep(paste0("Batch",1:B),nb)

start_time <- Sys.time()
snn.gr <- buildSNNGraph(sce, use.dimred="Corrected")
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

save.image(paste0("Comparison/",method, "_workspace.RData"))

#################################
# Apply Seurat to the LUAD Data #
#################################
rm(list=ls())
library(Seurat)
library(mclust) # For ARI
library(cluster) # For Silhouette
set.seed(12345)
method <- "Seurat"
proj <- "LUAD"
ver <- 1

##################
# Load LUAD Data #
##################
# Loading LUAD count data
countdata <- read.table(paste0("RawCountData/count_data_",proj,"_v",ver,".txt")) 

# Load dimension
dim <- unlist(read.table(paste0("RawCountData/dim_",proj,"_v",ver,".txt")))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

# Load gene_list
gene_list <- unlist(read.table(paste0("RawCountData/gene_list_",proj,"_v",ver,".txt"),stringsAsFactors = FALSE))

rownames(countdata) <- gene_list
colnames(countdata) <- rownames(metadata)

# Setting the Seurat Object
LUAD <- CreateSeuratObject(countdata, meta.data = metadata)
LUAD.list <- SplitObject(LUAD, split.by = "Protocol")


for (i in 1:length(LUAD.list)) {
  LUAD.list[[i]] <- NormalizeData(LUAD.list[[i]], verbose = FALSE)
  LUAD.list[[i]] <- FindVariableFeatures(LUAD.list[[i]], selection.method = "vst", nfeatures = G, 
                                         verbose = FALSE)
}

reference.list <- LUAD.list[c("10x", "CELseq2","Dropseq")]
LUAD.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = G)

LUAD.integrated <- IntegrateData(anchorset = LUAD.anchors, dims = 1:30)

DefaultAssay(LUAD.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
LUAD.integrated <- ScaleData(LUAD.integrated, verbose = FALSE)
LUAD.integrated <- RunPCA(LUAD.integrated, npcs = 30, verbose = FALSE)

LUAD.integrated <- FindNeighbors(LUAD.integrated, reduction = "pca", dims = 1:30, k.param = 20)
LUAD.integrated<-FindClusters(LUAD.integrated,resolution = 1.5)

w_Seurat <- LUAD.integrated$seurat_clusters
reorder <- match(rownames(metadata),names(w_Seurat))
w_Seurat<- w_Seurat[reorder]

Seurat_PCA <- LUAD.integrated@reductions$pca@cell.embeddings[reorder,]

Seurat_corrected <- NULL
for(g in 1:G){
  Seurat_corrected <- rbind(Seurat_corrected, LUAD.integrated@assays$integrated[g,])
}

#######
# ARI #
#######
ARI_Seurat <- adjustedRandIndex(metadata$CellType,w_Seurat)

save.image(paste0("Comparison/",method, "_workspace.RData"))

####################################
# Apply ZINB-WaVE to the LUAD Data #
####################################
rm(list=ls())
library(zinbwave)
library(mclust) # For ARI
library(cluster) # For Silhouette
set.seed(12345)
method <- "ZINBWaVE"
proj <- "LUAD"
ver <- 1

##################
# Load LUAD Data #
##################
# Loading LUAD count data
countdata <- read.table(paste0("RawCountData/count_data_",proj,"_v",ver,".txt")) 

# Load dimension
dim <- unlist(read.table(paste0("RawCountData/dim_",proj,"_v",ver,".txt")))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

# Load gene_list
gene_list <- unlist(read.table(paste0("RawCountData/gene_list_",proj,"_v",ver,".txt"),stringsAsFactors = FALSE))

rownames(countdata) <- gene_list
colnames(countdata) <- rownames(metadata)

##########################################
# Apply ZINB-WaVE to the Pancreatic Data #
##########################################
data_ZINBW <- as.matrix(countdata)

# Factorizing the batch indicators
batch_ind <- factor(rep(1:B,nb))

# Performing ZINB-WaVE
zinb_batch <- zinbFit(data_ZINBW, K = 10, X=model.matrix(~batch_ind), epsilon=1e3)

# Clustering
data_ZINBW <- SummarizedExperiment(assays = data_ZINBW)
assayNames(data_ZINBW)[1] <- "counts"
merged_zinb <- zinbwave(data_ZINBW, fitted_model = zinb_batch, K = 10, epsilon=1000)

library(Seurat)

seu <- as.Seurat(x = merged_zinb, counts = "counts", data = "counts")

seu
seu <- FindNeighbors(seu, reduction = "zinbwave",
                     dims = 1:10 #this should match K
)
seu <- FindClusters(object = seu)

w_ZINBWaVE <- seu$seurat_clusters

#######
# ARI #
#######
ARI_ZINBWaVE <- adjustedRandIndex(metadata$CellType,w_ZINBWaVE)

save.image(paste0("Comparison/",method, "_workspace.RData"))

########################################################
# Analyze the correction to the LUAD Data by Scanorama #
########################################################
rm(list=ls())
library(mclust) # For ARI
library(cluster) # For Silhouette
set.seed(12345)
method <- "scanorama"
proj <- "LUAD"
ver <- 1


##################
# Load LUAD Data #
##################
# load the dimentsion information
dim <- unlist(read.table(paste0("RawCountData/dim_",proj,"_v",ver,".txt")))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[1:B + 3]

# Loading the true cell type indicators
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

###################################
# load the inference by Scanorama #
###################################
scanorama_embedding <- read.table(paste0("Comparison/scanorama/scanorama_",proj,"_v",ver,"_integrated.txt"))
scanorama_corrected <- read.table(paste0("Comparison/scanorama/scanorama_",proj,"_v",ver,"_corrected.txt"))

clu_scanorama <- pam(scanorama_corrected, 3)
w_scanorama <- clu_scanorama$clustering

#######
# ARI #
#######
ARI_scanorama <- adjustedRandIndex(metadata$CellType, w_scanorama)

save.image(paste0("Comparison/",method, "_workspace.RData"))

###################################################
# Analyze the correction to the LUAD Data by scVI #
###################################################
rm(list=ls())
library(zinbwave)
library(mclust) # For ARI
library(cluster) # For Silhouette
set.seed(12345)
method <- "scVI"
proj <- "LUAD"
ver <- 1

##################
# Load LUAD Data #
##################
# load the dimentsion information
dim <- unlist(read.table(paste0("RawCountData/dim_",proj,"_v",ver,".txt")))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[1:B + 3]

# Loading the true cell type indicators
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

##############################
# load the inference by scVI #
##############################
w_scVI <- unlist(read.table(paste0("Comparison/scVI/scVI_",proj,"_v",ver,"_clusters.txt")))
scVI_corrected <- read.table(paste0("Comparison/scVI/scVI_",proj,"_v",ver,"_latent.txt"))

#######
# ARI #
#######
ARI_scVI <- adjustedRandIndex(metadata$CellType, w_scVI)

save.image(paste0("Comparison/",method, "_workspace.RData"))


##################
# Load workspace #
##################
load("BUSseq_workspace.RData")
comparison_list <- c("liger", "MNN", "scanorama", "scVI", "Seurat", "ZINBWaVE")
num_comparison <- length(comparison_list)

for(m in 1:num_comparison){
  load(paste0("Comparison/",comparison_list[m],"_workspace.RData"))
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

print(ARI_values)
