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
proj <- "LUAD"
ver <- 1
#######################
# Load Simulated Data #
#######################
# Working directory
setwd("DropletBased")

# Loading pancreatic count data
y_obs <- read.table(paste0("RawCountData/count_data_",proj,"_v",ver,".txt"))

# Load dimension
dim <- read.table(paste0("RawCountData/dim_",proj,"_v",ver,".txt"))
dim <- unlist(dim)
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3+1:B]

# Load metadata
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

# Load gene_list
gene_list <- unlist(read.table(paste0("RawCountData/gene_list_",proj,"_v",ver,".txt")))

if(!dir.exists("Image")){
  dir.create("Image")
}

###########################################
# Compare the BIC values for different Ks #
###########################################
k_sel <-2:6
N_k <- length(k_sel)

BIC.record <- rep(NA,N_k)
names(BIC.record) <- paste0("K=",k_sel)
for(k in k_sel){
  post_dir<-paste("Inference_K",k,sep="")
  if(dir.exists(post_dir)){
    # load BIC
    BIC<- unlist(read.table(paste(post_dir,"/BIC.txt",sep="")))
    BIC.record[k-k_sel[1] + 1] <- BIC[1]
    message(paste("Finish loading the results of ",proj,"_v",ver,"_K",k,"!/n",sep=""))
  }else{
    message(paste("The results of ",proj,"_v",ver,"_K",k," doesn't exist!/n",sep=""))
  }
}

jpeg(paste0("Image/",proj,"_v",ver,"_BIC.jpeg"),width = 800, height = 600)
par(mar=c(5.1,5.1,4.1,2.1))
plot(k_sel,BIC.record,xlab= "K",ylab = "BIC",type="n",cex.axis=2.5,cex.lab=3)
points(k_sel,BIC.record,type="b",pch=19,cex=3)
dev.off()

#####################################
# Load posterior inference of K = 5 #
#####################################
K <- 4

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
ARI_BUSseq <- adjustedRandIndex(metadata$CellType,w_BUSseq)

########################
# Set coloring to plot #
########################

#####set cell type colorings
color_by_celltype <- c(brewer.pal(4,"Set3"))
celltype_factor <- metadata$CellType

#####set batch colorings
color_by_batch <- c("#EB4334","#FBBD06","#35AA53")
batch_factor <- factor(metadata$Protocol)

# #####set cell type colorings
# color_by_celltype <- c(brewer.pal(3,"Set3"))
# celltype_factor <- metadata$CellType
# 
# #####set batch colorings
# color_by_batch <- c("#EB4334","#FBBD06","#35AA53")
# batch_factor <- factor(metadata$Protocol)

###############################
# Draw the Heatmap in Figure 3#
###############################
# color keys in the heatmap
# for mu.syn
color_key_mu_low <-  colorRampPalette(c("#0571b0","#F2F2F2"))
color_key_mu_high <-  colorRampPalette(c("#F2F2F2","#bd0026"))
color_key_mu <- c(color_key_mu_low(6),color_key_mu_high(6)[-1])

# # for nu.syn
# colorsChoice_nu_low <- colorRampPalette(c("#0571b0","#F2F2F2"))
# colorsChoice_nu_high <- colorRampPalette(c("#F2F2F2","#bd0026"))
# color_key_nu <- c(colorsChoice_nu_low(6),colorsChoice_nu_high(6)[-1])
# 
# for log_x log_y
colorsChoice_count <- colorRampPalette(c("#F2F2F2","black"))(11)

# The splitting points for log-scale mean expression levels and batch effects
break_mu <- seq(-1.5, 1.5, length.out = 12)
# break_nu <- seq(-4.2, 4.2, length.out = 12)
break_logcount <- seq(0,9,length.out = 12)

# Batch color bar for count data matrix
color_by_batch_count <- rep(color_by_batch,nb)

# Cell type color bar for count data matrix
color_by_celltype_count <- rep(NA,N)
celltype_est <- unlist(w.est) + 1
color_by_celltype_count <- color_by_celltype[celltype_est]


# Upper color bar for batch, lower color bar for cell type
col_annotation<-cbind(color_by_batch_count,color_by_celltype_count)
colnames(col_annotation) <- NULL

# reorder cells by cell type first and batch second
order_cells <-  order(w.est,rep(1:B,nb))
col_annotation <- col_annotation[order_cells,]

######################
# Plot the color key #
scale01 <- function(x, low = min(x), high = max(x)) {
  x <- (x - low)/(high - low)
  x
}

# color key for mu
z <- seq(min(break_mu), max(break_mu), length = length(color_key_mu))
jpeg(paste0("Image/heatmap_",proj,"_col_key_mu.jpg"),width = 360, height = 270, quality =100)
image(z = matrix(z, ncol = 1), col = color_key_mu, breaks = break_mu, xaxt = "n", yaxt = "n")
lv <- pretty(break_mu)
xv <- scale01(as.numeric(lv), min(break_mu), max(break_mu))
axis(1, at = xv, labels = lv, cex.axis = 3)
mtext(side = 1, "Value", line = 3.5, cex = 4)
title(main = "Color Key", cex.main = 4)
dev.off()

# color key for x_corrected
z <- seq(min(break_logcount), max(break_logcount), length = length(colorsChoice_count))
jpeg(paste0("Image/heatmap_",proj,"_col_key_logx.jpg"),width = 360, height = 270, quality =100)
image(z = matrix(z, ncol = 1), col = colorsChoice_count, breaks = break_logcount, xaxt = "n", yaxt = "n")
lv <- pretty(break_logcount)
xv <- scale01(as.numeric(lv), min(break_logcount), max(break_logcount))
axis(1, at = xv, labels = lv, cex.axis = 3)
mtext(side = 1, "Value", line = 3.5, cex = 4)
title(main = "Color Key", cex.main = 4)
dev.off()


############################################################################################
# The heatmap of the scaled estimated log-scale mean expression levels for intrinsic genes #
############################################################################################
log_mu_celltype_est <- alpha.est + beta.est
log_mu_intri <- log_mu_celltype_est[D.est == 1,]
# Sorting cell-type-specific expression levels 
# in the decreasing order of the proprotions of the first cell type
# because the true cell types are sorted in the same decreasing order
#order.est <- order(pi.est[1,],decreasing = T)


jpeg(paste0("Image/heatmap_",proj,"_est_log_mean_expression_levels.jpeg"),width = 1440, height = 720)
heatmap_mu_intri <- heatmap.3(log_mu_intri,
                    dendrogram = "none",#with cluster tree
                    Rowv = TRUE, Colv = FALSE,
                    labRow = FALSE, labCol = FALSE,
                    ColSideColors = color_by_celltype,
                    lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key 
                    #xlab = "Batch", ylab = "Intrinsic genes"
                    lhei=c(0.4,0.4,3.6),
                    col=color_key_mu,breaks = break_mu, key = FALSE)
dev.off()

###############################################################
# The heatmap of log-scale corrected read count data (Fig 3i) #
x_intri <- x_corrected[D.est==1,]
# Sharing the same double color bar and splitting points 
# as that of the underlying true count data
jpeg(paste0("Image/heatmap_",proj,"_log_corrected_count.jpeg"),width = 1440, height = 720)
heatmap.3(log1p(x_intri[rev(heatmap_mu_intri$rowInd),order_cells]),
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = col_annotation,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
          lhei=c(0.5,0.4,3.6),
          col=colorsChoice_count,breaks = break_logcount, key = FALSE)
dev.off()

##################################################################
# Draw t-SNE plots by batch and true cell types as Figure 4(c,d) #
##################################################################
set.seed(123)
all.dists.BUSseq <- as.matrix(dist(t(log1p(x_corrected[D.est==1,]))))
tsne_BUSseq_dist <- Rtsne(all.dists.BUSseq, is_distance=TRUE, perplexity = 30)

BUSseq_by_celltype<- paste0("Image/tsne_",proj,"_BUSseq_by_celltype.jpeg")
dat_frame <- data.frame(Var1 = tsne_BUSseq_dist$Y[,1], 
                        Var2 = tsne_BUSseq_dist$Y[,2], 
                        col = celltype_factor)

jpeg(BUSseq_by_celltype,width = 800, height = 400, quality = 100)
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

BUSseq_by_batch <- paste0("Image/tsne_",proj,"_BUSseq_by_batch.jpeg")
dat_frame <- data.frame(Var1 = tsne_BUSseq_dist$Y[,1], 
                        Var2 = tsne_BUSseq_dist$Y[,2], 
                        col = batch_factor)

jpeg(BUSseq_by_batch,width = 800, height = 400, quality = 100)
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

#####################
# Draw the PCA plot #
#####################
pca.BUSseq <- prcomp(t(log1p(x_corrected[D.est==1,])), rank=2)

BUSseq_PCA_by_celltype<- paste0("Image/pca_",proj,"_BUSseq_by_celltype.jpeg")
dat_frame <- data.frame(Var1 = pca.BUSseq$x[,1], 
                        Var2 = pca.BUSseq$x[,2], 
                        col = celltype_factor)

jpeg(BUSseq_PCA_by_celltype,width = 800, height = 400, quality = 100)
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

BUSseq_PCA_by_batch<- paste0("Image/pca_",proj,"_BUSseq_by_batch.jpeg")
dat_frame <- data.frame(Var1 = pca.BUSseq$x[,1], 
                        Var2 = pca.BUSseq$x[,2], 
                        col = batch_factor)

jpeg(BUSseq_PCA_by_batch,width = 800, height = 400, quality = 100)
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

# Store the workspace
save.image("BUSseq_workspace.RData")
save(ARI_BUSseq,tsne_BUSseq_dist,file = "BUSseq_results.RData")