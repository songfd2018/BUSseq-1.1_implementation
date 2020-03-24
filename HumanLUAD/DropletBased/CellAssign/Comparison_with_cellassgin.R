# Compare cellassign and BUSseq on the LUAD dataset
library(magrittr)
library(limma)
library(org.Hs.eg.db)
library(edgeR)
library(matrixStats)
library(pheatmap)
library(cellassign)
library(mclust)
library(MLmetrics)

setwd("DropletBased/CellAssign")

data(holik_data)
head(holik_data$counts[,1:2])
head(holik_data$cell_line)

entrez_map <- select(org.Hs.eg.db, 
                     as.character(rownames(holik_data$counts)), 
                     c("SYMBOL","ENSEMBL"), "ENTREZID")
gene_annotations <- entrez_map %>%
  dplyr::rename(GeneID=ENTREZID,
                Symbol=SYMBOL)

dge <- DGEList(counts = holik_data$counts, 
               group = holik_data$cell_line, 
               genes = gene_annotations, 
               remove.zeros = TRUE)
genes_to_keep <- rowSums(cpm(dge$counts) > 0.5) >= 2
dge_filt <- dge[genes_to_keep,]

dge_filt <- calcNormFactors(dge_filt, method="TMM")

# Differential expression

#We next perform differential expression using Limma Voom on a subset of 3 samples: HCC827, H2228, H1975:

dge_subset <- dge_filt[,dge_filt$samples$group %in% c("HCC827", "H2228", "H1975")]
design <- model.matrix(~ 0+dge_subset$samples$group)
colnames(design) <- levels(dge_subset$samples$group)
v <- voom(dge_subset, design)
fit <- lmFit(v, design)

# Next, fit constrasts to find differentially expressed genes between cell types:
contrast.matrix <- makeContrasts(H2228 - H1975, 
                                 HCC827 - H1975, 
                                 HCC827 - H2228, 
                                 levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Finally, compute gene summary statistics and filter to only significantly differentially expressed geens (FDR < 0.05):
  
tt <- topTable(fit2, n=Inf)
tt_sig <- tt %>%
  dplyr::filter(adj.P.Val < 0.05)

head(tt_sig)

# Marker gene derivation

#To derive marker genes, we first create a log fold change matrix using H1975 as the baseline expression:
  

lfc_table <- tt_sig[,c("H2228...H1975", "HCC827...H1975")]
lfc_table <- lfc_table %>%
  dplyr::mutate(H1975=0,
                H2228=H2228...H1975,
                HCC827=HCC827...H1975) %>%
  dplyr::select(H1975, H2228, HCC827)
rownames(lfc_table) <- tt_sig$GeneID

# Then, for each gene, we subtract the minimum log fold change, as 
# we care about overexpression of genes relative to some minimum 
# expression level, as this defines a marker gene:
  
lfc_table <- as.matrix(lfc_table)
lfc_table <- lfc_table - rowMins(lfc_table)
lfc_table <- as.data.frame(lfc_table)


# We now define a helper function for turning log fold changes into 
# a binary matrix. This takes a matrix and a threshold, and any values 
# less than or equal to the threshold are set to 0, and all others to 1:
  
binarize <- function(x, threshold) {
  x[x <= threshold] <- -Inf
  x[x > -Inf] <- 1
  x[x == -Inf] <- 0
  return(x)
}


# Find the biggest difference
maxdiffs <- apply(lfc_table, 1, function(x) max(diff(sort(x))))

#
thres_vals <- apply(lfc_table, 1, function(x) sort(x)[which.max(diff(sort(x)))])
expr_mat_thres <- plyr::rbind.fill(lapply(1:nrow(lfc_table), function(i) {
  binarize(lfc_table[i,], thres_vals[i])
}))
rownames(expr_mat_thres) <- rownames(lfc_table)
marker_gene_mat <- expr_mat_thres[(maxdiffs >= quantile(maxdiffs, c(.99))) 
                                  & (thres_vals <= log(2)),] %>%
  as.matrix


# Finally, we add back in gene symbols rather than entrez ids:
  

suppressMessages({
  symbols <- plyr::mapvalues(
    rownames(marker_gene_mat),
    from = gene_annotations$GeneID,
    to = gene_annotations$Symbol
  )
})

is_na <- is.na(symbols)

marker_gene_mat <- marker_gene_mat[!is_na,]
rownames(marker_gene_mat) <- symbols[!is_na]

#And there we have a marker gene matrix for our cell types:
  
head(marker_gene_mat)

pheatmap(t(marker_gene_mat))

# Turn the gene symbol to ENSEMBL ID
gene_index_match <- match(rownames(marker_gene_mat),entrez_map$SYMBOL)
rownames(marker_gene_mat) <- entrez_map$ENSEMBL[gene_index_match]

save(marker_gene_mat,file="LUAD_marker_genes.RData")

####################
# Apply CellAssign #
####################
library(SingleCellExperiment)
library(cellassign)
library(scran)

rm(list=ls())
set.seed(12345)
################
# load dataset #
################
load("../../sc_mixology-master/data/sincell_with_class.RData")
load("LUAD_marker_genes.RData")

# Obtain the common genes among three batches and calcualte the size factors
common_genes <- intersect(rownames(sce_sc_CELseq2_qc), intersect(rownames(sce_sc_10x_qc),rownames(sce_sc_Dropseq_qc)))

# raw count data of common genes
common_index <- match(common_genes, rownames(sce_sc_CELseq2_qc))
CountMat_CELseq2 <- counts(sce_sc_CELseq2_qc)[common_index,]
common_index <- match(common_genes, rownames(sce_sc_10x_qc))
CountMat_10x <- counts(sce_sc_10x_qc)[common_index,]
common_index <- match(common_genes, rownames(sce_sc_Dropseq_qc))
CountMat_Dropseq <- counts(sce_sc_Dropseq_qc)[common_index,]

# The true cell types are annotated for convenience
CellType <- c(colData(sce_sc_CELseq2_qc)$cell_line, colData(sce_sc_10x_qc)$cell_line, colData(sce_sc_Dropseq_qc)$cell_line)
Batch <- c(rep("CELseq2",ncol(sce_sc_CELseq2_qc)), rep("10x",ncol(sce_sc_10x_qc)), rep("Dropseq",ncol(sce_sc_Dropseq_qc)))

# construct the singlecellexperiment of three batches
sce_LUAD <- SingleCellExperiment(assays = list(counts = cbind(CountMat_CELseq2, CountMat_10x, CountMat_Dropseq)), 
                                 colData = list(Group = CellType, Batch = Batch))

# cell-specific size factor
sce_LUAD <- computeSumFactors(sce_LUAD)
sizefactors <- sizeFactors(sce_LUAD)

# Match the genes of each protocol to the marker gene
marker_gene_match <- match(rownames(marker_gene_mat),rownames(sce_LUAD))
marker_gene_index <- marker_gene_match[which(!is.na(marker_gene_match))]
marker_gene_mat <- marker_gene_mat[which(!is.na(marker_gene_match)),]

# batch covariates
batch_matrix <- model.matrix(~ 0 + (colData(sce_LUAD)$Batch=="10x") + (colData(sce_LUAD)$Batch=="Dropseq"))

# conduct cell assign
fit <- cellassign(exprs_obj = sce_LUAD[marker_gene_index,], 
                  marker_gene_info = marker_gene_mat, 
                  s = sizefactors,
                  X = batch_matrix,
                  learning_rate = 1e-2, 
                  shrinkage = TRUE,
                  verbose = FALSE)


print(fit)

# load metadata
metadata <- read.table("../RawCountData/metadata_LUAD_v1.txt")

# load the cell type indicators inferred by BUSseq
w_BUSseq <- unlist(read.table("../Inference_K3/w_est.txt"))
table(w_BUSseq,metadata$CellType)
celltype_BUSseq <- w_BUSseq
# Name each cell type
celltype_BUSseq[which(w_BUSseq==0)] <- "HCC827"
celltype_BUSseq[which(w_BUSseq==1)] <- "H2228"
celltype_BUSseq[which(w_BUSseq==2)] <- "H1975"

F1score <- rep(0,2)
Accu <- rep(0,2)
ARI <- rep(0,2)

F1score[1] <- F1_Score(metadata$CellType,celltype_BUSseq)
F1score[2] <- F1_Score(metadata$CellType,fit$cell_type)

Accu[1] <- Accuracy(metadata$CellType,celltype_BUSseq)
Accu[2] <- Accuracy(metadata$CellType,fit$cell_type)

ARI[1] <- adjustedRandIndex(metadata$CellType,celltype_BUSseq)
ARI[2] <- adjustedRandIndex(metadata$CellType,fit$cell_type)

results <- data.frame(F1 = F1score,
                      Acc = Accu,
                      ARI = ARI)
print(results)