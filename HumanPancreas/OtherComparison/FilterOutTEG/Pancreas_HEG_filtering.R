# Referring to the source code https://github.com/MarioniLab/MNN2017/blob/master/Pancreas

# This script prepares data for the pancreas analysis.
# It involves four publicly available datasets.

##########################################
##########################################
rm(list=ls())

setwd("OtherComparison/FilterOutTEG")

# we directly work on the normalized data in the "Data" folder
############
# GSE81076 #
############
# do a bit of tidying to make sure expression matrices are all conformable and
# only contain data with cell type assignments

# explicitly read in the study data
datah1 <- read.table("Data/GSE81076_SFnorm.tsv",
                     h=TRUE, sep="\t", stringsAsFactors=FALSE)

read1 <- read.table("Data/GSE81076_readcount.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)
# gene IDs are in the last column
genes1 <- datah1$gene_id
rownames(datah1) <- genes1
datah1 <- datah1[, 1:(dim(datah1)[2]-1)]

meta1 <- read.table("Data/GSE81076_marker_metadata.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)

# select meta data for which there is expression data
meta1 <- meta1[meta1$Sample %in% colnames(datah1), ]

# standardize the cell labels
## NB there are 112 samples with unassigned cell types, remove these
celltypes1 <- meta1$CellType
no.label1 <- meta1$Sample[meta1$CellType == ""]
datah1 <- datah1[, !colnames(datah1) %in% no.label1]
read1 <- read1[,!colnames(read1) %in% no.label1]
meta1 <- meta1[!meta1$Sample %in% no.label1, ]
celltypes1 <- celltypes1[celltypes1 != ""]
samples1 <- colnames(read1)

# check all dimensions match up
if(dim(datah1)[2] == dim(meta1)[1]) {dim(datah1)[2] == length(celltypes1)}


############
# GSE85241 #
############
# do a bit of tidying to make sure expression matrices are all conformable and
# only contain data with cell type assignments

# explicitly read in the study data
datah2 <- read.table("Data/GSE85241_SFnorm.tsv",
                     h=TRUE, sep="\t", stringsAsFactors=FALSE)

read2 <- read.table("Data/GSE85241_readcount.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)
# gene IDs are in the last column
genes2 <- datah2$gene_id
rownames(datah2) <- genes2
datah2 <- datah2[, 1:(dim(datah2)[2]-1)]

meta2 <- read.table("Data/GSE85241_marker_metadata.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)

# select meta data for which there is expression data
meta2 <- meta2[meta2$Sample %in% colnames(datah2), ]

# standardize the cell labels
## NB there are not any samples with unassigned cell types
celltypes2 <- meta2$CellType
no.label2 <- meta2$Sample[meta2$CellType == ""]
datah2 <- datah2[, !colnames(datah2) %in% no.label2]
read2 <- read2[,!colnames(read2) %in% no.label2]
meta2 <- meta2[!meta2$Sample %in% no.label2, ]
celltypes2 <- celltypes2[celltypes2 != ""]
samples2 <- colnames(read2)


# check all dimensions match up
if(dim(datah2)[2] == dim(meta2)[1]) {dim(datah2)[2] == length(celltypes2)}


############
# GSE86473 #
############
# do a bit of tidying to make sure expression matrices are all conformable and
# only contain data with cell type assignments

# explicitly read in the study data
datah3 <- read.table("Data/GSE86473_SFnorm.tsv",
                     h=TRUE, sep="\t", stringsAsFactors=FALSE)
read3 <- read.table("Data/GSE86473_readcount.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)
# gene IDs are in the last column
genes3 <- datah3$gene_id
rownames(datah3) <- genes3
datah3 <- datah3[, 1:(dim(datah3)[2]-1)]

meta3 <- read.table("Data/GSE86473_metadata.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)

# select meta data for which there is expression data
meta3 <- meta3[meta3$Sample %in% colnames(datah3), ]

# standardize the cell labels
## NB there are not any samples with unassigned cell types
celltypes3 <- meta3$CellType
no.label3 <- meta3$Sample[meta3$CellType == ""]
datah3 <- datah3[, !colnames(datah3) %in% no.label3]
read3 <- read3[,!colnames(read3) %in% no.label3]
meta3 <- meta3[!meta3$Sample %in% no.label3, ]
celltypes3 <- celltypes3[celltypes3 != ""]
samples3 <- colnames(read3)


###############
# E-MTAB-5061 #
###############
# do a bit of tidying to make sure expression matrices are all conformable and
# only contain data with cell type assignments
# explicitly read in the study data
datah4 <- read.table("Data/E-MTAB-5061_SFnorm.tsv",
                     h=TRUE, sep="\t", stringsAsFactors=FALSE)

read4 <- read.table("Data/E-MTAB-5061_readcount.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)
# gene IDs are in the last column
genes4 <- datah4$gene_id
rownames(datah4) <- genes4
datah4 <- datah4[, 1:(dim(datah4)[2]-1)]

meta4 <- read.table("Data/E-MTAB-5061_metadata.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)

# select meta data for which there is expression data
meta4 <- meta4[meta4$Sample %in% colnames(datah4), ]

# standardize the cell labels
## NB there are not any samples with unassigned cell types
celltypes4 <- meta4$CellType
no.label4 <- meta4$Sample[meta4$CellType == ""]
datah4 <- datah4[, !colnames(datah4) %in% no.label4]
read4 <- read4[,!colnames(read4) %in% no.label4]
meta4 <- meta4[!meta4$Sample %in% no.label4, ]
celltypes4 <- celltypes4[celltypes4 != ""]
samples4 <- colnames(read4)

# check all dimensions match up
if(dim(datah4)[2] == dim(meta4)[1]) {dim(datah4)[2] == length(celltypes4)}


# create one big meta data frame
all.meta <- do.call(rbind.data.frame, list("b1"=meta1[, c("Sample", "CellType", "Protocol", "Study")],
                                           "b2"=meta2[, c("Sample", "CellType", "Protocol", "Study")],
                                           "b3"=meta3[, c("Sample", "CellType", "Protocol", "Study")],
                                           "b4"=meta4[, c("Sample", "CellType", "Protocol", "Study")]))

# Take a look at the number of genes after qunality control and normalization
length(genes1)
length(genes2)
length(genes3)
length(genes4)

# Take logarithm and calculate the mean log-level
log_datah1 <- log1p(datah1)
mean_expr_datah1 <- apply(log_datah1,1,mean)

log_datah2 <- log1p(datah2)
mean_expr_datah2 <- apply(log_datah2,1,mean)

log_datah3 <- log1p(datah3)
mean_expr_datah3 <- apply(log_datah3,1,mean)

log_datah4 <- log1p(datah4)
mean_expr_datah4 <- apply(log_datah4,1,mean)

 
# Attain the most highly expressed genes
proportion <- 0.5
HE_ind_datah1 <- order(mean_expr_datah1) < proportion * length(mean_expr_datah1)
gene_sel_datah1 <- rownames(datah1)[which(HE_ind_datah1)]

HE_ind_datah2 <- order(mean_expr_datah2) < proportion * length(mean_expr_datah2)
gene_sel_datah2 <- rownames(datah2)[which(HE_ind_datah2)]

HE_ind_datah3 <- order(mean_expr_datah3) < proportion * length(mean_expr_datah3)
gene_sel_datah3 <- rownames(datah3)[which(HE_ind_datah3)]

HE_ind_datah4 <- order(mean_expr_datah4) < proportion * length(mean_expr_datah4)
gene_sel_datah4 <- rownames(datah4)[which(HE_ind_datah4)]

# Merge all four data matrices together based on a common set of gene IDs
common.genes <- intersect(gene_sel_datah1, intersect(gene_sel_datah2, intersect(gene_sel_datah3, gene_sel_datah4)))

mh1 <- match(common.genes, genes1)
mh2 <- match(common.genes, genes2)
mh3 <- match(common.genes, genes3)
mh4 <- match(common.genes, genes4)

raw.read.heg <- cbind(read1[mh1,], read2[mh2,], read3[mh3,], read4[mh4,])

# assign small weird cell types from GSE85241 and E=MTAB-5061 to 'other'
all.meta$CellType[grepl(all.meta$CellType, pattern="PP")] <- "Gamma"
all.meta$CellType[grepl(all.meta$CellType, pattern="Mesenchyme")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Co-ex")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Endo")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Epsi")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Mast")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="MHC")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Uncl")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Not")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="PSC")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="None")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Stellate")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="None")] <- "other"

nb <- c(ncol(read1), ncol(read2), ncol(read3), ncol(read4))

#store count data matrix
if(!dir.exists("RawCountData")){
  dir.create("RawCountData")
}


write.table(raw.read.heg,file="RawCountData/count_data_pancreas_v3.txt",row.names = FALSE, col.names = FALSE)

# txt version raw count data
mathed_index <- match(colnames(raw.read.heg),all.meta$Sample)
metadata <- data.frame(Protocol = all.meta$Protocol[mathed_index], Study = all.meta$Study[mathed_index], CellType = all.meta$CellType[mathed_index])
rownames(metadata) <- all.meta$Sample[mathed_index]

write.table(metadata,file="RawCountData/metadata_pancreas_v3.txt")

N <- ncol(raw.read.heg)
G <- nrow(raw.read.heg)
B <- length(nb)
write.table(c(N,G,B,nb),file="RawCountData/dim_pancreas_v3.txt",row.names = FALSE, col.names = FALSE)

write.table(rownames(raw.read.heg), file = "RawCountData/gene_list_pancreas_v3.txt",row.names = FALSE, col.names = FALSE)