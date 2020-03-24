#Apply BUSseq to the hematopoietic study.
rm(list=ls())
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)
library(gplots)
library(reshape2)

# coloring
library(scales)
require(WGCNA)
library(RColorBrewer)

library(biomaRt)
library(edgeR)

# plot slingshot
library(slingshot)
library(rgl)

library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

set.seed(123)
###################
# Load hemat Data #
###################
# Working directory
setwd("/your/working/directory/BUSseq_implementation-1.0/MouseHematopoietic/BUSseq/")

# Loading hemat count data
y_obs <- read.table("../RawCountData/count_data_hemat_v1.txt")

# Load dimension
dim <- read.table("../RawCountData/dim_hemat_v1.txt")
dim <- unlist(dim)
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3+1:B]

# Load metadata
metadata <- read.table("../RawCountData/metadata_hemat_v1.txt")

# Load gene_list
gene_list <- unlist(read.table("../RawCountData/gene_list_hemat_v1.txt",stringsAsFactors = F))

##########################
#load posterior inference#
##########################
K <- 6

# load w_est
w.est <- read.table("Inference_K6/w_est.txt")
w_BUSseq <- unlist(w.est)

# load alpha_est
# alpha.post <- as.matrix(read.big.matrix("alpha_post.txt",sep=" ",skip=num.burntin,type="double"))
alpha.est <- read.table("Inference_K6/alpha_est.txt")
alpha.est <- unlist(alpha.est)
# load beta_est
beta.est <- read.table("Inference_K6/beta_est.txt")
beta.est <- matrix(unlist(beta.est),G,K)
logmu.est<-beta.est+alpha.est

# load nu_est
nu.est <- read.table("Inference_K6/nu_est.txt")
nu.est <- matrix(unlist(nu.est),G,B)

# load delta_est
delta.est <- read.table("Inference_K6/delta_est.txt")
delta.est <- unlist(delta.est)
# plot(delta.est,col=rep(1:B,nb) + 1)


# load gamma_est
gamma.est <- read.table("Inference_K6/gamma_est.txt")
gamma.est <- matrix(unlist(gamma.est),B,2)


# load phi_est
phi.est <- read.table("Inference_K6/phi_est.txt")
phi.est <- matrix(unlist(phi.est),G,B)

# load pi_est
pi.est <- read.table("Inference_K6/pi_est.txt")
pi.est <- matrix(unlist(pi.est),B,K)
order.est<-order(pi.est[1,],decreasing = T)

# load p_est
p.est <- read.table("Inference_K6/p_est.txt")

# load tau0_est
tau0.est <- read.table("Inference_K6/tau0_est.txt")

# load PPI_est
PPI.est <- read.table("Inference_K6/PPI_est.txt")
D.est <- unlist(read.table("Inference_K6/IG_est.txt"))

# load X_imputed
x_imputed <- read.table("Inference_K6/x_imputed.txt")

############################
# Batch effects correction #
############################
adjusted_values = function(Truereads, .B, .nb, .N, .G, .K,
                          .alpha,.beta,.nu,.delta,.phi,.w){
    .w = .w + 1
    .b = unlist(sapply(seq_len(.B), function(b) rep(b, .nb[b])))
    #uniform random numbers
    global_u = array(runif(.N*.G),c(.G,.N))
    get_corrected_read = function(i){
        b = .b[i]
        # The following are vectors of length .G
        p_x <- pnbinom(Truereads[,i], size = .phi[,b], mu = exp(.alpha + .beta[,.w[i]] + .nu[,b] + .delta[i]))
        p_xminus1 <- pnbinom(Truereads[,i] - 1, size = .phi[,b], mu = exp(.alpha + .beta[,.w[i]] + .nu[,b] + .delta[i]))
        local_u = mapply(function(u) min(u, .999), global_u[,i]*(p_x - p_xminus1) + p_xminus1)
        return(qnbinom(local_u, size = .phi[,1], mu = exp(.alpha + .beta[,.w[i]])))
    }
    return(mapply(get_corrected_read, seq_len(.N)))
    
}

start_time<-Sys.time()
print("Calculate corrected read counts:")
x_corrected<-adjusted_values(x_imputed, B, nb, N, G, K,
                             alpha.est,beta.est,nu.est,delta.est,phi.est,w_BUSseq)
write.table(x_corrected, file = "x_corrected.txt",row.names = FALSE, col.names = FALSE)
x_intrinsic <- x_corrected[D.est==1,]

end_time<-Sys.time()
running_time<-end_time-start_time
print(running_time)

####################
# Pathway analysis #
####################
intri_gene <- gene_list[D.est==1]
write.table(intri_gene, file = paste0("intrinsic_gene_list.txt"), row.names = FALSE,col.names = FALSE,quote = FALSE)

#######
# ARI #
#######
ARI_BUSseq <- adjustedRandIndex(metadata$CellType,w_BUSseq)

kclust_BUSseq <- pam(t(log1p(x_intrinsic)), k = 7)
w_BUSseq_kmeans <- kclust_BUSseq$clustering
ARI_BUSseq_kmeans <- adjustedRandIndex(w_BUSseq_kmeans,metadata$CellType)

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

# cell type est color bar
color_by_celltype_est<-c("#1A9850", # 5
                         "#feb24c", # 3
                         "#41b6c4", # 0
                         "#fd8d3c", # 1
                         "#D7301F", # 4
                         "#FFFF99" # 2
                     )

celltype_est_factor <- factor(w_BUSseq + 1, levels = c("6","4","1","2","5","3"))

##############################################################
# Draw t-SNE plots by batch and true cell types Figure 5(a,b)#
##############################################################
set.seed(123)
all.dists.BUSseq <- as.matrix(dist(t(log1p(x_corrected[D.est==1,]))))
tsne_BUSseq_dist <- Rtsne(all.dists.BUSseq, is_distance=TRUE, perplexity = 30)

BUSseq_by_celltype<- "Image/tsne_hemat_BUSseq_by_celltype.jpg"
dat_frame <- data.frame(Var1 = tsne_BUSseq_dist$Y[,1], 
                        Var2 = tsne_BUSseq_dist$Y[,2], 
                        col = celltype_factor)

jpeg(BUSseq_by_celltype,width = 800, height = 600, quality = 100)
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

BUSseq_by_batch <- "Image/tsne_hemat_BUSseq_by_batch.jpg"
dat_frame <- data.frame(Var1 = tsne_BUSseq_dist$Y[,1], 
                        Var2 = tsne_BUSseq_dist$Y[,2], 
                        col = batch_factor)

jpeg(BUSseq_by_batch,width = 800, height = 600, quality = 100)
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

################################
# Draw the PCA plot Figure 5(c)#
################################
pca.BUSseq <- prcomp(t(log1p(x_corrected[D.est==1,])), rank=2)

BUSseq_PCA<- "Image/pca_hemat_BUSseq.jpeg"
#plot_by_celltype(BUSseq_PCA, pca.BUSseq$x, xlab = "PC 1", ylab = "PC 2")
dat_frame <- data.frame(Var1 = pca.BUSseq$x[,1], 
                        Var2 = pca.BUSseq$x[,2], 
                        col = celltype_factor)

jpeg(BUSseq_PCA,width = 800, height = 600, quality = 100)
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

######################################################
# Draw PCA plots by estimated cell types Figure 6(b) #
######################################################
pca_frame <- data.frame(PC1 = pca.BUSseq$x[,1], PC2 = pca.BUSseq$x[,2],
                        Celltype = celltype_est_factor)
jpeg("Image/pca_hemat_BUSseq_est.jpeg",width = 800, height = 800, quality = 100)
p <- ggplot(pca_frame, aes(x=PC1, y=PC2, colour=Celltype)) +
  geom_point(size = 4) + theme_classic() +
  guides(colour = guide_legend(override.aes = list(size = 12))) +
  scale_colour_manual(values = alpha(color_by_celltype_est,0.6), labels = paste0("Cluster ",1:K)) +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

###############################################
# Draw the heatmap of marker genes Figure 6(e)#
###############################################
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

key_genes <- unlist(read.table("hemat_key_genes.txt",stringsAsFactors = FALSE, quote = ""))
key_gene_symbol <- getBM(attributes=c('ensembl_gene_id',
                                      'external_gene_name'),
                         filters = 'mgi_symbol',
                         values = key_genes,
                         mart = ensembl)

# adjust the log-scale mean expression levels to be non-negative
mu_adjust <- log(1+exp(alpha.est + beta.est))
rownames(mu_adjust) <- gene_list
mu_adjust <- mu_adjust[,c(6,4,1,2,5,3)]
colnames(mu_adjust) <- paste0("Cluster ",1:K)
mu_adjust_scaled <- t(scale(t(mu_adjust)))

matched_key_genes <- match(key_gene_symbol$ensembl_gene_id,gene_list)
mu_key_genes <- mu_adjust_scaled[matched_key_genes[!is.na(matched_key_genes)],]
rownames(mu_key_genes) <- key_gene_symbol$external_gene_name[which(!is.na(matched_key_genes))]

# plot by ggplot2
hc.cols <- hclust(dist(mu_key_genes))
gene_order <- hc.cols$order
#mu.frame <- melt(mu_key_genes[gene_order,])
mu.frame <- melt(mu_key_genes[c(5, 34, 19, 26, 13,
                                6, 17, 4, 21, 35,
                                22, 33, 10, 37, 20, 
                                36,  28,  1, 15, 29, 
                                16, 14,  2, 23,  3,  
                                8, 18, 24, 25, 27,
                                32, 9, 30, 31, 11,  
                                7, 12),])
colnames(mu.frame) <- c("Gene","CellType","Mu")
jpeg("Image/hemat_heatmap_key_genes_in_cell_paper.jpg",width = 700, height = 800, quality = 100)
p <- ggplot(mu.frame, aes(CellType, Gene)) + geom_tile(aes(fill = Mu),colour = "white") + 
  scale_fill_gradient2(low = "#0571b0",
                        midpoint = 0,
                        mid = "#fcfcfc",
                        high = "#bd0026",
                        space="Lab") +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Cell type by BUSseq") + labs(fill = "Scaled\nExpression\nLevels") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 32), #angle = 45),
        axis.text.y = element_text(face = "bold.italic", #color = "blue", 
                                   size = 18),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=24,face="bold"),
        legend.title = element_text(size = 36,face="bold"))
p
dev.off()

#############################################################################################################
# Calculate the correlation between the estimated cell type effects and mean expression levels in GSE116177 #
#############################################################################################################
Haemopedia <- read.table(file = "../RawCountData/GSE116177_Haemopedia-Mouse-RNASeq_raw.txt",header = TRUE,stringsAsFactors = FALSE)
rownames(Haemopedia) <- Haemopedia$geneId
Haemopedia <- Haemopedia[,-1]

# normlization
TMM_factor <- calcNormFactors(Haemopedia, method = "TMM")
N_Haemopedia <- ncol(Haemopedia)
Haemopedia_normlized <- Haemopedia
log_Haemopedia_normlized <- Haemopedia
for(i in 1:N_Haemopedia){
  Haemopedia_normlized[,i] <- Haemopedia[,i] / TMM_factor[i]
}
log_Haemopedia_normlized <- log1p(Haemopedia_normlized)

# match the shared genes between BUSseq and Haemopedia genes
gene_matched_Haemopedia <- match(gene_list,rownames(log_Haemopedia_normlized))
gene_matched_BUSseq <- which(!is.na(gene_matched_Haemopedia))
gene_matched_Haemopedia <- gene_matched_Haemopedia[gene_matched_BUSseq]

mu.est <- alpha.est[gene_matched_BUSseq] + beta.est[gene_matched_BUSseq,]
rownames(mu.est) <- gene_list[gene_matched_BUSseq]
log_mu.est <- log1p(exp(mu.est))

# adjust the order of cell type

log_mu.est <- log_mu.est[,c(6,4,1,2,5,3)]
colnames(log_mu.est) <- paste0("Cluster ",1:K)
celltype_est_adjusted <- w_BUSseq
celltype_est_adjusted[which(w_BUSseq==0)] <- 3
celltype_est_adjusted[which(w_BUSseq==1)] <- 4
celltype_est_adjusted[which(w_BUSseq==2)] <- 6
celltype_est_adjusted[which(w_BUSseq==3)] <- 2
celltype_est_adjusted[which(w_BUSseq==4)] <- 5
celltype_est_adjusted[which(w_BUSseq==5)] <- 1

log_Haemopedia_hvgs <- log_Haemopedia_normlized[gene_matched_Haemopedia,]
log_Haemopedia_hvgs <- as.matrix(log_Haemopedia_hvgs)
# assign names for all cells
colnames(log_Haemopedia_hvgs) <- c("MegTPO-1", "MegTPO-2", "MegSCF-1", "MegSCF-2", "BasoCult-1",	"PreCFUE-1",
                                   "CFUE-1", "CFUE-2", "PreCFUE-2", "LSK-1", "NeutBM-1", "NeutBM-2", "Eo-1",
                                   "CD8T-1", "CD8T-2", "CD4T-1", "CD4T-2", "Mac-1", "Mac-2", "LSK-2", "BasoCult-2",
                                   "EoP-1", "EoP-2", "GMP-1", "GMP-2", "Mast-1", "MonoPB-1", "MonoPB-2",
                                   "CLP-1", "Mast-2", "BCell-1", "BCell-2", "MonoBM-1", "MonoBM-2", "MEP-1", "MEP-2",
                                   "CMP-1", "CMP-2", "InfMono-1", "InfMono-2", "MPP-1", "MPP-2", "STHSC-1", "STHSC-2",
                                   "EryBlPO-1", "EryBlPB-1", "EryBlPB-2", "NeutPB-1", "Eo-2", "Eo-3", "Retic-1",
                                   "CLP-2", "NeutPB-2", "Retic-2", "Retic-3", "BasoSpl-1", "EryBlPO-2", "EryBlPO-3", "EryBlPB-3",
                                   "BasoBM-1", "BasoBM-2", "MemCD4T-1", "EffCD4T-1", "EffCD4T-2", "NveCD4T-1", "NveCD4T-2", 
                                   "RegT-1", "NK-1", "NK-2", "MemCD8T-1", "NveCd8T-1", "NveCd8T-2", "NeutPB-3", "NeutPB-4",
                                   "RegT-2", "MemCD4T-2", "MegCD8T-1", "pDC-1", "pDC-2", "cDC2-1", "cDC2-2", "GMPSiglecF-1", 
                                   "EoSSCHi-1", "EoSSCLo-1", "EoSSCLo-2", "GMPSiglecF-2", "EoSSCHi-2", "EoSSCLo-3",
                                   "EoSSCHi-3", "GMPSiglecF-3", "EoSSCHi-4", "EoSSCHi-5", "EoCult-1", "EoCult-2", 
                                   "MacCult-1", "MacCult-2", "MacCult-3", "GMPSiglecF-4", "EoSSCHi-6",
                                   "EoSSCLo-4", "EoSSCLo-5", "GMPSiglecF-5")

# take average of expression levels for the selected cell types
log_Haemopedia_hvgs_selected <- NULL
celltype_list <- NULL
# LSK
celltype_name <- "LSK"
select_index <- grep(celltype_name,colnames(log_Haemopedia_hvgs))
length(select_index)
express_mean <- apply(log_Haemopedia_hvgs[,select_index],1,mean)
log_Haemopedia_hvgs_selected <- cbind(log_Haemopedia_hvgs_selected, express_mean)
celltype_list <- c(celltype_list, celltype_name)
# STHSC
celltype_name <- "STHSC"
select_index <- grep(celltype_name,colnames(log_Haemopedia_hvgs))
length(select_index)
express_mean <- apply(log_Haemopedia_hvgs[,select_index],1,mean)
log_Haemopedia_hvgs_selected <- cbind(log_Haemopedia_hvgs_selected, express_mean)
celltype_list <- c(celltype_list, celltype_name)
# MPP
celltype_name <- "MPP"
select_index <- grep(celltype_name,colnames(log_Haemopedia_hvgs))
length(select_index)
express_mean <- apply(log_Haemopedia_hvgs[,select_index],1,mean)
log_Haemopedia_hvgs_selected <- cbind(log_Haemopedia_hvgs_selected, express_mean)
celltype_list <- c(celltype_list, celltype_name)
# CLP
celltype_name <- "CLP"
select_index <- grep(celltype_name,colnames(log_Haemopedia_hvgs))
length(select_index)
express_mean <- apply(log_Haemopedia_hvgs[,select_index],1,mean)
log_Haemopedia_hvgs_selected <- cbind(log_Haemopedia_hvgs_selected, express_mean)
celltype_list <- c(celltype_list, celltype_name)
# CMP
celltype_name <- "CMP"
select_index <- grep(celltype_name,colnames(log_Haemopedia_hvgs))
length(select_index)
express_mean <- apply(log_Haemopedia_hvgs[,select_index],1,mean)
log_Haemopedia_hvgs_selected <- cbind(log_Haemopedia_hvgs_selected, express_mean)
celltype_list <- c(celltype_list, celltype_name)
# MEP
celltype_name <- "MEP"
select_index <- grep(celltype_name,colnames(log_Haemopedia_hvgs))
length(select_index)
express_mean <- apply(log_Haemopedia_hvgs[,select_index],1,mean)
log_Haemopedia_hvgs_selected <- cbind(log_Haemopedia_hvgs_selected, express_mean)
celltype_list <- c(celltype_list, celltype_name)
# GMP
celltype_name <- "GMP-"
select_index <- grep(celltype_name,colnames(log_Haemopedia_hvgs))
length(select_index)
express_mean <- apply(log_Haemopedia_hvgs[,select_index],1,mean)
log_Haemopedia_hvgs_selected <- cbind(log_Haemopedia_hvgs_selected, express_mean)
celltype_name <- "GMP"
celltype_list <- c(celltype_list, celltype_name)
colnames(log_Haemopedia_hvgs_selected) <- celltype_list

mu_key <- log_mu.est[rownames(mu.est) %in% key_gene_symbol$ensembl_gene_id,]
Haemopedia_key <- log_Haemopedia_hvgs_selected[rownames(mu.est) %in% key_gene_symbol$ensembl_gene_id,]

mu_key_scaled <- t(scale(t(mu_key)))
Haemopedia_key_scaled <- t(scale(t(Haemopedia_key)))
cor_BUSseq_Haemopedia_key_scaled <- cor(mu_key_scaled, Haemopedia_key_scaled)
cor_intri_key_scaled <- melt(cor_BUSseq_Haemopedia_key_scaled)
colnames(cor_intri_key_scaled) <- c("CellType_BUSseq","CellType_Haemopedia_selected","Correlation")
jpeg("Image/Normlized_expression_correlation_key_genes_selected_celltype_after_scaling.jpg",width = 700, height = 800, quality = 100)
p <- ggplot(cor_intri_key_scaled, aes(x = CellType_BUSseq, y = CellType_Haemopedia_selected, fill = Correlation)) + 
  geom_tile(aes(),colour = "white") +
  scale_fill_gradient2(low = "#41ab5d",
                       midpoint = 0,
                       mid = "#fcfcfc",
                       high = "#fc4e2a",
                       space="Lab") +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Cell type by BUSseq") + ylab("Cell type in the Haemopedia") +
  theme(axis.text.x = element_text(face = "bold",
                                   size = 32), 
        axis.text.y = element_text(face = "bold",
                                   size = 32),
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=24,face="bold"),
        legend.title = element_text(size = 32,face="bold"))
p
dev.off()


##############################################
# Draw the expression levels of marker genes #
##############################################
# plot heatmap of marker genes
marker_genes <- c("Apoe","Gata2","Ctse","Ebf1","Vpreb1", "Igll1")
marker_genes_symbol <- getBM(attributes=c('ensembl_gene_id',
                                          'external_gene_name'),
                             filters = 'mgi_symbol',
                             values = marker_genes,
                             mart = ensembl)

# Apoe gene
ensembl_id <- marker_genes_symbol$ensembl_gene_id[which(marker_genes_symbol$external_gene_name == "Apoe")]
Apoe_index <- which(gene_list==ensembl_id)

# Gata2 gene
ensembl_id <- marker_genes_symbol$ensembl_gene_id[which(marker_genes_symbol$external_gene_name == "Gata2")]
Gata2_index <- which(gene_list==ensembl_id)

# Ctse gene for MEP
ensembl_id <- marker_genes_symbol$ensembl_gene_id[which(marker_genes_symbol$external_gene_name == "Ctse")]
Ctse_index <- which(gene_list==ensembl_id)

# Igll1 gene for B-biased CLP
ensembl_id <- marker_genes_symbol$ensembl_gene_id[which(marker_genes_symbol$external_gene_name == "Igll1")]
Igll1_index <- which(gene_list==ensembl_id)


marker_genes_frame <- data.frame(Cell_Index = 1:N, Protocol = metadata$Protocol, 
                                 PCA_intri_x =  pca.BUSseq$x[,1], PCA_intri_y =  pca.BUSseq$x[,2],
                                 Study = metadata$Study, Celltype = factor(celltype_est_adjusted),
                                 Celltype_est = w_BUSseq,
                                 Apoe = log1p(unlist(x_corrected[Apoe_index,])),
                                 Gata2 = log1p(unlist(x_corrected[Gata2_index,])), 
                                 Ctse = log1p(unlist(x_corrected[Ctse_index,])),
                                 Igll1 = log1p(unlist(x_corrected[Igll1_index,])))

# Expression levels of APOE gene indicating different stages of erythropoiesis
jpeg(filename = "Image/Hemat_APOE_from_CMP_to_MEP.jpeg",width = 700, height = 700,quality = 100)
ggplot(marker_genes_frame, aes(x=PCA_intri_x, y=PCA_intri_y, colour=Apoe)) +
  geom_text(aes(label=Celltype), size = 8) + theme_classic() +
  scale_color_gradient2(low = "#0571b0",
                       midpoint = 0,
                       mid = "#fcfcfc",
                       high = "#bd0026",
                       space="Lab") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold",
                                   size = 36),
        axis.text.y = element_text(face = "bold",
                                   size = 36),
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=24,face="bold"),
        legend.title = element_text(size = 36,face="bold.italic"),
        panel.background = element_rect(colour = "black",size = 2))
dev.off()

# Expression levels of GATA2 gene indicating different stages of erythropoiesis
jpeg(filename = "Image/Hemat_GATA2_from_CMP_to_MEP.jpeg",width = 700, height = 700,quality = 100)
ggplot(marker_genes_frame, aes(x=PCA_intri_x, y=PCA_intri_y, colour=Gata2,
                               shape=Celltype)) +
  geom_text(aes(label=Celltype), size = 8) + theme_classic() +
  scale_color_gradient2(low = "#0571b0",
                        midpoint = 0,
                        mid = "#fcfcfc",
                        high = "#bd0026",
                        space="Lab") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", 
                                   size = 36),
        axis.text.y = element_text(face = "bold",
                                   size = 36),
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=24,face="bold"),
        legend.title = element_text(size = 36,face="bold.italic"),
        panel.background = element_rect(colour = "black",size = 2))
dev.off()

# Higher expression levels of CTSE gene in cluster 5
jpeg(filename = "Image/Hemat_CTSE_for_MEP.jpeg",width = 700, height = 700,quality = 100)
ggplot(marker_genes_frame, aes(x=PCA_intri_x, y=PCA_intri_y, colour=Ctse,
                               shape=Celltype)) +
  geom_text(aes(label=Celltype), size = 8) + theme_classic() +
  scale_color_gradient2(low = "#0571b0",
                        midpoint = 0,
                        mid = "#fcfcfc",
                        high = "#bd0026",
                        space="Lab") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", 
                                   size = 36),
        axis.text.y = element_text(face = "bold",
                                   size = 36),
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=24,face="bold"),
        legend.title = element_text(size = 36,face="bold.italic"),
        panel.background = element_rect(colour = "black",size = 2))
dev.off()

# Higher expression levels of Igll1 gene in cluster 3
jpeg(filename = "Image/Hemat_Igll1_for_LMPP.jpeg",width = 700, height = 700,quality = 100)
ggplot(marker_genes_frame, aes(x=PCA_intri_x, y=PCA_intri_y, colour=Igll1,
                               shape=Celltype)) +
  geom_text(aes(label=Celltype), size = 8) + theme_classic() +
  scale_color_gradient2(low = "#0571b0",
                        midpoint = 0,
                        mid = "#fcfcfc",
                        high = "#bd0026",
                        space="Lab") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", 
                                   size = 36),
        axis.text.y = element_text(face = "bold", 
                                   size = 36),
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=24,face="bold"),
        legend.title = element_text(size = 36,face="bold.italic"),
        panel.background = element_rect(colour = "black",size = 2))
dev.off()

##################################################
# box plot of cell type effects vs batch effects #
##################################################
beta_nu_est <- cbind(beta.est,nu.est)
colnames(beta_nu_est) <- c(paste0("CellType",1:K),paste0("Batch",1:B))

# exclude the first cell type and first batch
beta_nu_est <- beta_nu_est[,-c(1,K+1)]

beta_nu_est_frame <- melt(beta_nu_est)
colnames(beta_nu_est_frame) <- c("Gene","Var","Value")

jpeg(filename = "Image/hemat_boxplot_of_beta_and_nu.jpg",width = 900, height = 600,quality = 100)
ggplot(beta_nu_est_frame, aes(x = Var, y = Value, fill = Var)) + 
  geom_boxplot() + theme_bw() +
  geom_vline(xintercept = 5.5, linetype="dashed", size = 1) +
  ylim(-5,5)+
  scale_fill_manual(values = c(color_by_celltype_est[-1],color_by_batch[-1])) +
  theme(axis.text.x = element_text(face = "bold", 
                                   size = 36, angle = 90),
        axis.text.y = element_text(face = "bold", 
                                   size = 36),
        axis.title=element_text(size=44,face="bold"),
        legend.position = "none",
        axis.title.x=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
dev.off()

#####################################################
# Slingshot of hematopoietic development trajectory #
#####################################################
pca.BUSseq.3d <- prcomp(t(log1p(x_intrinsic)), rank=3)
celltype_est_factor <- factor(w_BUSseq + 1, levels = c("6","4","1","2","5","3"))
BUSseq_linage_3d <- slingshot(pca.BUSseq.3d$x, clusterLabels = celltype_est_factor, start.clus = "6", end.clus = c("1","3","5"))
plot3d(pca.BUSseq.3d$x, col = color_by_celltype_est[celltype_est_factor])
plot3d(BUSseq_linage_3d, lwd = 3, add = TRUE)

# Store the workspace
save(ARI_BUSseq, ARI_BUSseq_kmeans, tsne_BUSseq_dist,x_intrinsic, file = "BUSseq_results.RData")