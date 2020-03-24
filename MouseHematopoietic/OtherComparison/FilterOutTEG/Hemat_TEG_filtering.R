# Referring to the source code https://github.com/MarioniLab/MNN2017/blob/master/Haematopoiesis/prepareData.R

# This script prepares data for the hematopoiesis analysis.
# It involves two publicly available datasets.

##########################################
##########################################
rm(list=ls())

setwd("OtherComparison/FilterOutTEG")

if(!dir.exists("RawData")){
  dir.create("RawData")
}

###################
###################
## prepareData.R ##
###################
###################

# Download and read the counts, metadata of Nestorowa et al. 2016
fname <- "RawData/GSE81682_HTSeq_counts.txt.gz"
if (!file.exists(fname)) { download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81682&format=file&file=GSE81682%5FHTSeq%5Fcounts%2Etxt%2Egz", fname) }
dataF <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
dataF <- as.matrix(dataF)
dim(dataF)

fname <- "RawData/metaF.txt"
if (!file.exists(fname)) { download.file("http://blood.stemcells.cam.ac.uk/data/all_cell_types.txt", fname) }
metaF <- read.table(fname, stringsAsFactors = FALSE, header=TRUE, check.names=FALSE)
metainds <- match(colnames(dataF), rownames(metaF))
missing.meta <- is.na(metainds)
metaF <- metaF[metainds,] # This will contain NA's... which is okay, at this point, to preserve length.

# Defining the cell type based on the metadata.
metatypeF <- rep("other", nrow(metaF))
for (col in rev(colnames(metaF))) { # reverse, so earlier columns end up overwriting later ones.
  chosen <- metaF[,col]==1
  metatypeF[chosen] <- sub("[0-9]?_.*", "", col)
}
metatypeF[metatypeF=="ESLAM"] <- "HSPC"

# Filling in metadata from the cell sorting label, if metadata was missing.
metatypeF[missing.meta] <- sub("_.*", "", colnames(dataF)[missing.meta])
metatypeF[metatypeF=="LT-HSC"] <- "LTHSC"
metatypeF[metatypeF=="Prog"] <- "other"
colnames(dataF)<-metatypeF

# Perform size factor normalization within this data set.
library(scran)
high.abF <- scater::calcAverage(dataF) > 1
clustF <- quickCluster(dataF, method="igraph", subset.row=high.abF)
sizeF <- computeSumFactors(dataF, cluster=clustF, subset.row=high.abF)
dataF2 <- t(t(dataF)/sizeF)

# Cleaning up memory.
gc() 

##########################################
##########################################

# Download and read the counts and meta data of Paul et al. 2015
fname <- "RawData/GSE72857_umitab.txt.gz"
if (!file.exists(fname)) { download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE72857&format=file&file=GSE72857%5Fumitab%2Etxt%2Egz", fname) }
dataA <- read.table(fname, header=TRUE, row.names=1)
metaA <- read.csv2("RawData/MAP.csv",sep=",",stringsAsFactors = FALSE, head=TRUE, row.names=1)
dim(dataA)

# Only selecting cells that are in the metadata.
metainds <- match(rownames(metaA), colnames(dataA))
dataA <- dataA[,metainds]
dataA <- as.matrix(dataA)

# Organizing cell type labels.
metatypeA <- character(nrow(metaA))
metatypeA[metaA[,1]<7] <- "ERY"
metatypeA[metaA[,1]>6 & metaA[,1]<12] <- "CMP"
metatypeA[metaA[,1]>11] <- "GMP"
colnames(dataA) <- metatypeA

# Perform size factor normalization within this data set.
high.abA <- scater::calcAverage(dataA) > 1
clustA <- quickCluster(dataA, method="igraph", subset.row=high.abA)
sizeA <- computeSumFactors(dataA, cluster=clustA, subset.row=high.abA)
dataA2 <- t(t(dataA)/sizeA)

# Cleaning up memory.
gc() 

save.image("RawCountDataAfterNormalization.RData")

################################################################
################################################################
## Select the 10% most highly expressed gene in each datasets ##
################################################################
################################################################
# Take logarithm and calculate the mean log-level
log_dataA2 <- log1p(dataA2)
mean_expr_dataA <- apply(log_dataA2,1,mean)

log_dataF2 <- log1p(dataF2)
mean_expr_dataF <- apply(log_dataF2,1,mean)

# Attain the most highly expressed genes
proportion <- 0.5
HE_ind_dataA <- order(mean_expr_dataA) < proportion * length(mean_expr_dataA)
gene_sel_dataA <- rownames(dataA2)[which(HE_ind_dataA)]

HE_ind_dataF <- order(mean_expr_dataF) < proportion * length(mean_expr_dataF)
gene_sel_dataF <- rownames(dataF2)[which(HE_ind_dataF)]


# Pull down IDs from BioMaRt.
library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org" )
out <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), values = gene_sel_dataF, mart = mart,filters = "ensembl_gene_id")

match_genes <- pmatch(out$mgi_symbol, gene_sel_dataA) # partial, due to use of concatenated gene symbols.
length(which(!is.na(match_genes)))

ind_out <- which(!is.na(match_genes))
interest_gene <- cbind(out$ensembl_gene_id[ind_out],gene_sel_dataA[match_genes[ind_out]])

# Selecting features that are TEGs and present in both data sets.
mF <- match(interest_gene[,1], rownames(dataF2))
mA <- match(interest_gene[,2],rownames(dataA2))

dataA3 <- dataA[mA,]
dataF3 <- dataF[mF,]
rownames(dataA3) <- rownames(dataF3)

save.image(file="RawCountData_hemat_v3.Rdata")

if(!dir.exists("RawCountData")){
  dir.create("RawCountData")
}

nb <- c(ncol(dataA3),ncol(dataF3))

colnames(dataA3)[which(colnames(dataA3)=="ERY")] <- "MEP"

metadata <- data.frame(Protocol = rep(c("MARS-seq","Smart-seq2"),nb), 
                       Study = rep(c("GSE72857","GSE81682"),nb),
                       CellType = c(colnames(dataA3),colnames(dataF3)))

write.table(metadata,file="RawCountData/metadata_hemat_v3.txt")

hemat_count <- cbind(dataA3,dataF3)
write.table(hemat_count,file="RawCountData/count_data_hemat_v3.txt",row.names = FALSE, col.names = FALSE)

N <- ncol(hemat_count)
G <- nrow(hemat_count)
B <- length(nb)
write.table(c(N,G,B,nb),file="RawCountData/dim_hemat_v3.txt",row.names = FALSE, col.names = FALSE)

write.table(rownames(dataA3), file ="RawCountData/gene_list_hemat_v3.txt",row.names = FALSE, col.names = FALSE)