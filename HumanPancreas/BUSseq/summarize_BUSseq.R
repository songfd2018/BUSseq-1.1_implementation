#Apply BUSseq to the pancreatic study.
rm(list=ls())
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)
library(reshape2) # For melt
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

# coloring
library(scales)
library(WGCNA)
library(RColorBrewer)

set.seed(12345)

########################
# Load Pancreatic Data #
########################
# Working directory
setwd("BUSseq")

# Loading pancreatic count data
y_obs <- read.table("../RawCountData/count_data_pancreas_v1.txt")

# Load dimension
dim <- read.table("../RawCountData/dim_pancreas_v1.txt")
dim <- unlist(dim)
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3+1:B]

# Load metadata
metadata <- read.table("../RawCountData/metadata_pancreas_v1.txt")

# Load gene_list
gene_list <- unlist(read.table("../RawCountData/gene_list_pancreas_v1.txt",stringsAsFactors = F))

##########################
#load posterior inference#
##########################
K <- 8

# load w_est
w.est <- read.table("../Inference_K8/w_est.txt")
w_BUSseq <- unlist(w.est)

# load alpha_est
# alpha.post <- as.matrix(read.big.matrix("alpha_post.txt",sep=" ",skip=num.burntin,type="double"))
alpha.est <- read.table("../Inference_K8/alpha_est.txt")
alpha.est <- unlist(alpha.est)
# load beta_est
beta.est <- read.table("../Inference_K8/beta_est.txt")
beta.est <- matrix(unlist(beta.est),G,K)
logmu.est<-beta.est+alpha.est

# load nu_est
nu.est <- read.table("../Inference_K8/nu_est.txt")
nu.est <- matrix(unlist(nu.est),G,B)

# load delta_est
delta.est <- read.table("../Inference_K8/delta_est.txt")
delta.est <- unlist(delta.est)
plot(delta.est,col=rep(1:B,nb) + 1)


# load gamma_est
gamma.est <- read.table("../Inference_K8/gamma_est.txt")
gamma.est <- matrix(unlist(gamma.est),B,2)


# load phi_est
phi.est <- read.table("Inference_K8/phi_est.txt")
phi.est <- matrix(unlist(phi.est),G,B)

# load pi_est
pi.est <- read.table("Inference_K8/pi_est.txt")
pi.est <- matrix(unlist(pi.est),B,K)
order.est<-order(pi.est[1,],decreasing = T)

# load p_est
p.est <- read.table("Inference_K8/p_est.txt")

# load tau0_est
tau0.est <- read.table("Inference_K8/tau0_est.txt")

# load PPI_est
PPI.est <- read.table("Inference_K8/PPI_est.txt")
D.est <- unlist(read.table("Inference_K8/IG_est.txt"))

# load X_imputed
x_imputed <- read.table("Inference_K8/x_imputed.txt")


####################
# Pathway analysis #
####################
intri_gene <- gene_list[D.est==1]
write.table(intri_gene, file = paste0("intrinsic_gene_list.txt"), row.names = FALSE,col.names = FALSE,quote = FALSE)

############################
# Batch effects correction #
############################
adjusted_values = function(Truereads, .B, .nb, .N, .G, .K,
                          .alpha,.beta,.nu,.delta,.phi,.w){
    .w = .w + 1
    .b = unlist(sapply(seq_len(.B), function(b) rep(b, .nb[b])))
    #uniform random numbers are obtained outside of parallel code to handle seed issue.
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
write.table(x_corrected, file = "x_corrected.txt", row.names = FALSE, col.names = FALSE)
x_intrinsic <- x_corrected[D.est==1,]
end_time<-Sys.time()
running_time<-end_time-start_time
print(running_time)



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
color_by_celltype <- c(brewer.pal(6,"Set3"),"black")
color_by_celltype[c(1,6)] <- color_by_celltype[c(6,1)]
celltype_factor <- factor(metadata$CellType, levels = c("Alpha", "Beta", "Gamma", "Delta", "Acinar", "Ductal", "other"))

color_by_celltype_est <- c(brewer.pal(7,"Set3"),"black")
color_by_celltype_est <- color_by_celltype_est[c(4,1,6,2,3,8,5,7)]
celltype_est_factor <- factor(w_BUSseq)

#####set batch colorings
color_by_batch <- c("#EB4334","#FBBD06","#35AA53","#4586F3")
batch_factor <- factor(rep(1:B,nb))

#################################################
# Draw t-SNE plots by batch and true cell types #
#################################################
set.seed(123)
all.dists.BUSseq <- as.matrix(dist(t(log1p(x_corrected[D.est==1,]))))
tsne_BUSseq_dist <- Rtsne(all.dists.BUSseq, is_distance=TRUE, perplexity = 30)

BUSseq_by_celltype<- "Image/tsne_pancreas_BUSseq_by_celltype.jpg"
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

BUSseq_by_batch <- "Image/tsne_pancreas_BUSseq_by_batch.jpg"
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

BUSseq_by_celltype_est<- "Image/tsne_pancreas_BUSseq_by_celltype_est.jpg"
dat_frame <- data.frame(Var1 = tsne_BUSseq_dist$Y[,1], 
                        Var2 = tsne_BUSseq_dist$Y[,2], 
                        col = celltype_est_factor)
jpeg(BUSseq_by_celltype,width = 800, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_celltype_est,0.6)) +
  xlab("tSNE 1") + ylab("tSNE 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 44), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 44),#, angle = 45))
        axis.title=element_text(size=48,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.position = "none")
dev.off()

##################
# Draw PCA plots #
##################
pca.BUSseq <- prcomp(t(log1p(x_corrected[D.est==1,])), rank=2)
BUSseq_PCA<- "Image/pca_pancreas_BUSseq.jpeg"
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

BUSseq_PCA_est<- "Image/pca_pancreas_BUSseq_est.jpeg"
dat_frame <- data.frame(Var1 = pca.BUSseq$x[,1], 
                        Var2 = pca.BUSseq$x[,2], 
                        col = celltype_est_factor)

jpeg(BUSseq_PCA,width = 800, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_celltype_est,0.6)) +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.position = "none")
dev.off()

##############################################
# Draw the expression levels of marker genes #
##############################################
# GCG for alpha cell
GCG_index <-  which(gene_list=="GCG")

# INS for beta cells
INS_index <-  which(gene_list=="INS")

# SST for delta cells
SST_index <-  which(gene_list=="SST")

# PPY for gamma cells
PPY_index <-  which(gene_list=="PPY")

# PRSS1 for acinar cells
PRSS1_index <-  which(gene_list=="PRSS1")

# KRT19 for ductal cells
KRT19_index <-  which(gene_list=="KRT19")

marker_genes <- data.frame(Cell_Index = 1:N, Protocol = metadata$Protocol, 
                           tSNE1 =  tsne_BUSseq_dist$Y[,1], tSNE2 =  tsne_BUSseq_dist$Y[,2],
                           Study = metadata$Study, Celltype = metadata$CellType, 
                           GCG = log1p(unlist(x_corrected[GCG_index,])),
                           INS = log1p(unlist(x_corrected[INS_index,])), 
                           SST = log1p(unlist(x_corrected[SST_index,])),
                           PPY = log1p(unlist(x_corrected[PPY_index,])),
                           PRSS1 = log1p(unlist(x_corrected[PRSS1_index,])),
                           KRT19 = log1p(unlist(x_corrected[KRT19_index,])))

# Higher expression levels of GCG gene in alpha cells
jpeg(filename = "Image/Pancreas_GCG_for_alpha.jpeg",width = 800, height = 600,quality = 100)
ggplot(marker_genes, aes(x=tSNE1, y=tSNE2, colour=GCG,
                         shape=(Celltype=="Alpha"))) +
  geom_point(size=3) + theme_classic() +
  scale_color_gradient2(low = "#0571b0",
                        midpoint = 0,
                        mid = "#ffffcc",
                        high = "#bd0026",
                        space="Lab") +
  scale_shape_manual(values = c(4, 16), breaks = c(TRUE,FALSE)) +
  guides(shape = guide_legend(order = 1, 
                              override.aes = list(size = 12))) + # change symbol size in the legend 
  labs(shape = "Is Alpha Cell") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
    size = 36), #angle = 45),
    axis.text.y = element_text(face = "bold", #color = "blue", 
      size = 36),#, angle = 45))
    axis.title=element_text(size=44,face="bold"),
    legend.text=element_text(size=28,face="bold"),
    legend.title = element_text(size = 36,face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
dev.off()

# Higher expression levels of INS gene in beta cells
jpeg(filename = "Image/Pancreas_INS_for_beta.jpeg",width = 800, height = 600,quality = 100)
ggplot(marker_genes, aes(x=tSNE1, y=tSNE2, colour=INS,
                         shape=(Celltype=="Beta"))) +
  geom_point(size=3) + theme_classic() +
  scale_color_gradient2(low = "#0571b0",
                        midpoint = 0,
                        mid = "#ffffcc",
                        high = "#bd0026",
                        space="Lab") +
  scale_shape_manual(values = c(4, 16), breaks = c(TRUE,FALSE)) +
  guides(shape = guide_legend(order = 1, 
                              override.aes = list(size = 12))) + # change symbol size in the legend 
  labs(shape = "Is Beta Cell") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"),
        panel.background = element_rect(colour = "black",size = 2))
dev.off()

# Higher expression levels of SST gene in delta cells
jpeg(filename = "Image/Pancreas_SST_for_delta.jpeg",width = 800, height = 600,quality = 100)
ggplot(marker_genes, aes(x=tSNE1, y=tSNE2, colour=SST,
                         shape=(Celltype=="Delta"))) +
  geom_point(size=3) + theme_classic() +
  scale_color_gradient2(low = "#0571b0",
                        midpoint = 0,
                        mid = "#ffffcc",
                        high = "#bd0026",
                        space="Lab") +
  scale_shape_manual(values = c(4, 16), breaks = c(TRUE,FALSE)) +
  guides(shape = guide_legend(order = 1, 
                              override.aes = list(size = 12))) + # change symbol size in the legend 
  labs(shape = "Is Gamma Cell") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"),
        panel.background = element_rect(colour = "black",size = 2))
dev.off()

# Higher expression levels of PPY gene in gamma cells
jpeg(filename = "Image/Pancreas_PPY_for_gamma.jpeg",width = 800, height = 600,quality = 100)
ggplot(marker_genes, aes(x=tSNE1, y=tSNE2, colour=PPY,
                         shape=(Celltype=="Gamma"))) +
  geom_point(size=3) + theme_classic() +
  scale_color_gradient2(low = "#0571b0",
                        midpoint = 0,
                        mid = "#ffffcc",
                        high = "#bd0026",
                        space="Lab") +
  scale_shape_manual(values = c(4, 16), breaks = c(TRUE,FALSE)) +
  guides(shape = guide_legend(order = 1, 
                              override.aes = list(size = 12))) + # change symbol size in the legend 
  labs(shape = "Is Delta Cell") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"),
        panel.background = element_rect(colour = "black",size = 2))
dev.off()

# Higher expression levels of PRSS1 gene in acinar cells
jpeg(filename = "Image/Pancreas_PRSS1_for_acinar.jpeg",width = 800, height = 600,quality = 100)
ggplot(marker_genes, aes(x=tSNE1, y=tSNE2, colour=PRSS1,
                         shape=(Celltype=="Acinar"))) +
  geom_point(size=3) + theme_classic() +
  scale_color_gradient2(low = "#0571b0",
                        midpoint = 0,
                        mid = "#ffffcc",
                        high = "#bd0026",
                        space="Lab") +
  scale_shape_manual(values = c(4, 16), breaks = c(TRUE,FALSE)) +
  guides(shape = guide_legend(order = 1, 
                              override.aes = list(size = 12))) + # change symbol size in the legend 
  labs(shape = "Is Acinar Cell") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"),
        panel.background = element_rect(colour = "black",size = 2))
dev.off()

# Higher expression levels of KRT19 gene in ductal cells
jpeg(filename = "Image/Pancreas_KRT19_for_ductal.jpeg",width = 800, height = 600,quality = 100)
ggplot(marker_genes, aes(x=tSNE1, y=tSNE2, colour=KRT19,
                         shape=(Celltype=="Ductal"))) +
  geom_point(size=3) + theme_classic() +
  scale_color_gradient2(low = "#0571b0",
                        midpoint = 0,
                        mid = "#ffffcc",
                        high = "#bd0026",
                        space="Lab") +
  scale_shape_manual(values = c(4, 16), breaks = c(TRUE,FALSE)) +
  guides(shape = guide_legend(order = 1, 
                              override.aes = list(size = 12))) + # change symbol size in the legend 
  labs(shape = "Is Ductal Cell") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"),
        panel.background = element_rect(colour = "black",size = 2))
dev.off()

##################################################
# box plot of cell type effects vs batch effects #
##################################################
beta_nu_est <- cbind(beta.est,nu.est)
colnames(beta_nu_est) <- c(paste0("CellType ",1:K),paste0("Batch ",1:B))

# exclude the first cell type and first batch
beta_nu_est <- beta_nu_est[,-c(1,K+1)]

beta_nu_est_frame <- melt(beta_nu_est)
colnames(beta_nu_est_frame) <- c("Gene","Var","Value")

jpeg(filename = "Image/pancreas_boxplot_of_beta_and_nu.jpg",width = 900, height = 600,quality = 100)
ggplot(beta_nu_est_frame, aes(x = Var, y = Value, fill = Var)) + 
  geom_boxplot() + theme_bw() +
  geom_vline(xintercept = 7.5, linetype="dashed", size = 1) +
  ylim(-5,5)+
  scale_fill_manual(values = c(color_by_celltype_est[-1],color_by_batch[-1])) +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36, angle = 90),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        legend.position = "none",
        axis.title.x=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
dev.off()

#plot heatmap
# cluster genes
# mu_nu_est <- cbind(alpha.est + beta.est,nu.est)
# colnames(mu_nu_est) <- c(paste0("CellType",1:K),paste0("Batch",1:B))
# mu_nu_est <- mu_nu_est[,-(K+1)]

ord <- hclust( dist(beta_nu_est, method = "euclidean"), method = "ward.D" )$order
beta_nu_est_frame$Gene <- factor(beta_nu_est_frame$Gene, levels = ord)

jpeg("Image/heatmap_pancreas_beta_nu_comparison.jpg",width = 600, height = 900, quality = 100)
p <- ggplot(beta_nu_est_frame, aes(x = Var, y = Gene, fill = Value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#41ab5d",
                       midpoint = 0,
                       mid = "#fcfcfc",
                       high = "#fc4e2a",
                       space="Lab") +
  theme(axis.text.x = element_text(angle = 90), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) +
  #xlab("Cell type by BUSseq") + ylab("Cell type in the Haemopedia") +
  theme(axis.text.x = element_text(face = "bold",
                                   size = 32), 
        #axis.text.y = element_text(face = "bold",
        #                           size = 32),
        axis.title=element_text(size=44,face="bold"),
        legend.text=element_text(size=24,face="bold"),
        legend.title = element_text(size = 32,face="bold"))
p
dev.off()

##################################################################################
# Plot the heatmap of raw count data and corrected read count data for all genes #
##################################################################################
col_annotation <- rep(color_by_batch,nb)
colorsChoice_grey <- colorRampPalette(c("#EEEEEE","black"))
color_key <- colorsChoice_grey(5)
logy <- log1p(y_obs)

break_read <- seq(0,4,length.out = 6)
jpeg("Image/heatmap_pancreas_log_raw_count.jpeg",width = 1440, height = 1080)
heatmap.3(logy,
          dendrogram = "none",#with cluster tree
          Rowv = TRUE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
          lhei=c(0.6,0.4,3.6),
          col = color_key, breaks = break_read,
          ColSideColors = col_annotation, key = FALSE)
dev.off()

logx <- log1p(x_corrected)
jpeg("Image/heatmap_pancreas_log_correct_count.jpeg",width = 1440, height = 1080)
heatmap.3(logx,
          dendrogram = "none",#with cluster tree
          Rowv = TRUE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
          lhei=c(0.6,0.4,3.6),
          col = color_key, breaks = break_read,
          ColSideColors = col_annotation, key = FALSE)
dev.off()

######################
# Plot the color key #
scale01 <- function(x, low = min(x), high = max(x)) {
  x <- (x - low)/(high - low)
  x
}
# color key
z <- seq(min(break_read), max(break_read), length.out = 6)
jpeg("Image/heatmap_pancreas_col_key_y.jpg",width = 360, height = 270, quality =100)
image(z = matrix(z, ncol = 1), col = color_key, breaks = break_read, xaxt = "n", yaxt = "n")
lv <- pretty(break_read)
xv <- scale01(as.numeric(lv), min(break_read), max(break_read))
axis(1, at = xv, labels = lv, cex.axis = 3)
mtext(side = 1, "Value", line = 3.5, cex = 4)
title(main = "Color Key", cex.main = 4)
dev.off()


###########################################################################
# plot the heatmap of corrected read count matrix for all intrinsic genes #
###########################################################################
cell_est_factor <- factor(w_BUSseq, levels = c(1,7,0,3,4,6,2,5))
col_annotation <- cbind(rep(color_by_batch,nb),c(brewer.pal(7,"Set3"),"black")[cell_est_factor])
cell_order <- order(cell_est_factor)
col_annotation <- col_annotation[cell_order,]

colorsChoice_low <- colorRampPalette(c("blue","#F2F2F2"))
colorsChoice_high <- colorRampPalette(c("#F2F2F2","red"))
color_key <- c(colorsChoice_low(6),colorsChoice_high(6)[-1])

x_scaled <- log1p(x_corrected)[which(D.est==1),cell_order]
x_scaled <- t(scale(t(x_scaled)))
break_x_scaled <- seq(-2,2,length.out = 12)
jpeg("Image/heatmap_pancreas_log_corrected_count.jpeg",width = 1440, height = 1080)
heatmap.3(x_scaled,
          dendrogram = "none",#with cluster tree
          Rowv = TRUE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
          lhei=c(0.6,0.4,3.6),
          col = color_key, breaks = break_x_scaled,
          ColSideColors = col_annotation, key = FALSE)
dev.off()

######################
# Plot the color key #
scale01 <- function(x, low = min(x), high = max(x)) {
  x <- (x - low)/(high - low)
  x
}
# color key
z <- seq(min(break_x_scaled), max(break_x_scaled), length = length(color_key))
jpeg("Image/heatmap_pancreas_col_key_x_scaled.jpg",width = 360, height = 270, quality =100)
image(z = matrix(z, ncol = 1), col = color_key, breaks = break_x_scaled, xaxt = "n", yaxt = "n")
lv <- pretty(break_x_scaled)
xv <- scale01(as.numeric(lv), min(break_x_scaled), max(break_x_scaled))
axis(1, at = xv, labels = lv, cex.axis = 3)
mtext(side = 1, "Value", line = 3.5, cex = 4)
title(main = "Color Key", cex.main = 4)
dev.off()


# Store the workspace
save(ARI_BUSseq, ARI_BUSseq_kmeans, tsne_BUSseq_dist,x_intrinsic, file = "BUSseq_results.RData")
save.image("BUSseq_workspace.RData")