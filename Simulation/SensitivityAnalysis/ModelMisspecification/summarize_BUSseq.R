#Apply BUSseq to the simulation study.
rm(list=ls())
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)

# coloring
library(scales)
library(WGCNA)
library(RColorBrewer)

# heatmap.3 referring https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

set.seed(12345)
ver <- 2

#######################
# Load Simulated Data #
#######################
# Working directory
setwd("SensitivityAnalysis/ModelMisspecification")

# Loading pancreatic count data
y_obs <- read.table(paste0("RawCountData/count_data_simulation_v",ver,".txt"))

# Load dimension
dim <- read.table(paste0("RawCountData/dim_simulation_v",ver,".txt"))
dim <- unlist(dim)
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3+1:B]

# Load metadata
metadata <- read.table(paste0("RawCountData/metadata_simulation_v",ver,".txt"))

# Load gene_list
gene_list <- paste0("gene-",1:G)

#####################################
# Load posterior inference of K = 5 #
#####################################
K <- 5

# load w_est
w.est <- read.table(paste("Inference_K",K,"/w_est.txt",sep=""))
w_BUSseq <- unlist(w.est)

# load alpha_est
# alpha.post <- as.matrix(read.big.matrix("alpha_post.txt",sep=" ",skip=num.burntin,type="double"))
alpha.est <- read.table(paste("Inference_K",K,"/alpha_est.txt",sep=""))
alpha.est <- unlist(alpha.est)
# load beta_est
beta.est <- read.table(paste("Inference_K",K,"/beta_est.txt",sep=""))
beta.est <- matrix(unlist(beta.est),G,K)
logmu.est<-beta.est+alpha.est

# load nu_est
nu.est <- read.table(paste("Inference_K",K,"/nu_est.txt",sep=""))
nu.est <- matrix(unlist(nu.est),G,B)

# load delta_est
delta.est <- read.table(paste("Inference_K",K,"/delta_est.txt",sep=""))
delta.est <- unlist(delta.est)

# load gamma_est
gamma.est <- read.table(paste("Inference_K",K,"/gamma_est.txt",sep=""))
gamma.est <- matrix(unlist(gamma.est),B,2)


# load phi_est
phi.est <- read.table(paste("Inference_K",K,"/phi_est.txt",sep=""))
phi.est <- matrix(unlist(phi.est),G,B)

# load pi_est
pi.est <- read.table(paste("Inference_K",K,"/pi_est.txt",sep=""))
pi.est <- matrix(unlist(pi.est),B,K)
order.est<-order(pi.est[1,],decreasing = T)

# load p_est
p.est <- read.table(paste("Inference_K",K,"/p_est.txt",sep=""))

# load tau0_est
tau0.est <- read.table(paste("Inference_K",K,"/tau0_est.txt",sep=""))

# load PPI_est
PPI.est <- read.table(paste("Inference_K",K,"/PPI_est.txt",sep=""))
D.est <- unlist(read.table(paste("Inference_K",K,"/IG_est.txt",sep="")))

# load X_imputed
x_imputed <- read.table(paste("Inference_K",K,"/x_imputed.txt",sep=""))


################################################
# Prepare intrinsic genes for Pathway analysis #
################################################
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
end_time<-Sys.time()
running_time<-end_time-start_time
print(running_time)

#######
# ARI #
#######
ARI_BUSseq <- adjustedRandIndex(metadata$celltype,w_BUSseq)

########################
# Set coloring to plot #
########################

#set cell type colorings
# batch color bar
color_by_batch<-c("#EB4334","#FBBD06","#35AA53","#4586F3")
batch_factor <- factor(paste0("Batch ",rep(1:B,nb)))
batch.cols<-color_by_batch[rep(1:B,nb)]

# cell type color bar
color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978","#8097D3","#9C7ACE")
celltype_factor <- factor(paste0("Celltype ",metadata$celltype))

###############################
# Draw the Heatmap in Figure 3#
###############################
# color keys in the heatmap
# for mu.syn
color_key_mu <-  colorRampPalette(c("#F2F2F2","black"))(11)

# for nu.syn
colorsChoice_nu_low <- colorRampPalette(c("#0571b0","#F2F2F2"))
colorsChoice_nu_high <- colorRampPalette(c("#F2F2F2","#bd0026"))
color_key_nu <- c(colorsChoice_nu_low(6),colorsChoice_nu_high(6)[-1])

# for log_x log_y
colorsChoice_count <- colorRampPalette(c("#F2F2F2","black"))(11)

# The splitting points for log-scale mean expression levels and batch effects
break_mu <- seq(-0.5, 4, length.out = 12)
break_nu <- seq(-4.2, 4.2, length.out = 12)
break_logcount <- seq(0,9,length.out = 12)

##########################################################################
# The heatmap of the estimated log-scale mean expression levels (Fig 3f) #
log_mu_celltype_est <- alpha.est + beta.est

# Sorting cell-type-specific expression levels 
# in the decreasing order of the proprotions of the first cell type
# because the true cell types are sorted in the same decreasing order
order.est <- order(pi.est[1,],decreasing = T)

jpeg("Image/heatmap_simulation_est_log_mean_expression_levels.jpg",width = 1080, height = 1440)
heatmap.3(log_mu_celltype_est[, order.est],
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = color_by_celltype,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key 
          lhei=c(0.4,0.4,3.6),
          col=color_key_mu,breaks = break_mu, key = FALSE)
dev.off()

##################################################
# The heatmap of location batch effects (Fig 3g) #
batch_effects_est <- nu.est

jpeg("Image/heatmap_simulation_est_batch_effects.jpg",width = 1080, height = 1440)
heatmap.3(batch_effects_est,
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = color_by_batch,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key 
          lhei=c(0.6,0.4,3.6),
          col=color_key_nu,breaks = break_nu, key = FALSE)
dev.off()

#####################################################################
# The heatmap of log-scale underlying true read count data (Fig 3h) #
log_imputed_count_est <- log1p(x_imputed)

# Batch color bar for count data matrix
color_by_batch_count <- rep(color_by_batch,nb)

# Cell type color bar for count data matrix
color_by_celltype_count <- rep(NA,N)
celltype_est <- unlist(w.est) + 1
for(i in 1:N){
  color_by_celltype_count[i] <- color_by_celltype[which(celltype_est[i]==order.est)]
}


# Upper color bar for batch, lower color bar for cell type
col_annotation<-cbind(color_by_celltype_count,color_by_batch_count)
colnames(col_annotation) <- NULL

jpeg("Image/heatmap_simulation_log_imputed_count.jpg",width = 1080, height = 1440)
heatmap.3(log_imputed_count_est,
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = col_annotation,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
          lhei=c(0.6,0.4,3.6),
          col=colorsChoice_count,breaks = break_logcount, key = FALSE)
dev.off()


###############################################################
# The heatmap of log-scale corrected read count data (Fig 3i) #

# Sharing the same double color bar and splitting points 
# as that of the underlying true count data
jpeg("Image/heatmap_simulation_log_corrected_count.jpg",width = 1080, height = 1440)
heatmap.3(log1p(x_corrected),
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = col_annotation,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
          lhei=c(0.5,0.4,3.6),
          col=colorsChoice_count,breaks = break_logcount, key = FALSE)
dev.off()

####################################################################################################################
# The scatter plot of the estimated cell-specific size factors versus the true cell-specific size factors (Fig 3e) #

# Loading the true cell-specific size factors
cell_effect_true_vec <- unlist(read.table("True_para/delta_syn.txt"))

# Obtaining the estimated cell-specific size factors
cell_effect_est <- delta.est

jpeg("Image/scatter_simulation_cell_size_factors.jpg",width = 900, height = 600)
par(mar = c(5.1,6.1,4.1,2.1)) 
plot(cell_effect_true_vec,cell_effect_est,xlab = expression(paste("True ",delta)), ylab = expression(paste("Estimated ",delta)),type="n",ylim = c(-4,2),xlim = c(-4,2),cex.axis = 3, cex.lab = 3)
points(cell_effect_true_vec,cell_effect_est,pch = 1,cex=3)
abline(a=0,b=1,lty=3,cex=3)
dev.off()

##################################################################
# Draw t-SNE plots by batch and true cell types as Figure 4(c,d) #
##################################################################
set.seed(123)
all.dists.BUSseq <- as.matrix(dist(t(log1p(x_corrected[D.est==1,]))))
tsne_BUSseq_dist <- Rtsne(all.dists.BUSseq, is_distance=TRUE, perplexity = 30)

BUSseq_by_celltype<- "Image/tsne_simulation_BUSseq_by_celltype.jpg"
dat_frame <- data.frame(Var1 = tsne_BUSseq_dist$Y[,1], 
                        Var2 = tsne_BUSseq_dist$Y[,2], 
                        col = celltype_factor)

jpeg(BUSseq_by_celltype,width = 900, height = 600, quality = 100)
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

BUSseq_by_batch <- "Image/tsne_simulation_BUSseq_by_batch.jpg"
dat_frame <- data.frame(Var1 = tsne_BUSseq_dist$Y[,1], 
                        Var2 = tsne_BUSseq_dist$Y[,2], 
                        col = batch_factor)

jpeg(BUSseq_by_batch,width = 900, height = 600, quality = 100)
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

#####################
# Draw the PCA plot #
#####################
pca.BUSseq <- prcomp(t(log1p(x_corrected[D.est==1,])), rank=2)

BUSseq_PCA<- "Image/pca_simulation_BUSseq.jpg"
dat_frame <- data.frame(Var1 = pca.BUSseq$x[,1], 
                        Var2 = pca.BUSseq$x[,2], 
                        col = celltype_factor)

jpeg(BUSseq_PCA,width = 900, height = 600, quality = 100)
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
save.image("BUSseq_workspace.RData")
save(ARI_BUSseq,tsne_BUSseq_dist,file = "BUSseq_results.RData")