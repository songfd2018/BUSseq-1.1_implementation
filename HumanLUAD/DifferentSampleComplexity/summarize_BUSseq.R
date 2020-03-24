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

setwd("DifferentSampleComplexity")
set.seed(12345)
proj <- "LUAD"
for(ver in 2:3){

#######################
# Load Simulated Data #
#######################
if(ver == 2){
  folder <- "three_mix/"
}else{
  folder <- "five_mix/"
}
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
gene_list <- read.table(paste0("RawCountData/gene_list_",proj,"_v",ver,".txt"))

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
K <- 3

# load w_est
w.est <- read.table(paste(folder,"Inference_K",K,"/w_est.txt",sep=""))
w_BUSseq <- unlist(w.est)

# load alpha_est
# alpha.post <- as.matrix(read.big.matrix("alpha_post.txt",sep=" ",skip=num.burntin,type="double"))
alpha.est <- read.table(paste(folder,"Inference_K",K,"/alpha_est.txt",sep=""))
alpha.est <- unlist(alpha.est)
# load beta_est
beta.est <- read.table(paste(folder,"Inference_K",K,"/beta_est.txt",sep=""))
beta.est <- matrix(unlist(beta.est),G,K)
logmu.est<-beta.est+alpha.est

# load nu_est
nu.est <- read.table(paste(folder,"Inference_K",K,"/nu_est.txt",sep=""))
nu.est <- matrix(unlist(nu.est),G,B)

# load delta_est
delta.est <- read.table(paste(folder,"Inference_K",K,"/delta_est.txt",sep=""))
delta.est <- unlist(delta.est)

# load gamma_est
gamma.est <- read.table(paste(folder,"Inference_K",K,"/gamma_est.txt",sep=""))
gamma.est <- matrix(unlist(gamma.est),B,2)


# load phi_est
phi.est <- read.table(paste(folder,"Inference_K",K,"/phi_est.txt",sep=""))
phi.est <- matrix(unlist(phi.est),G,B)

# load pi_est
pi.est <- read.table(paste(folder,"Inference_K",K,"/pi_est.txt",sep=""))
pi.est <- matrix(unlist(pi.est),B,K)
order.est<-order(pi.est[1,],decreasing = T)

# load p_est
p.est <- read.table(paste(folder,"Inference_K",K,"/p_est.txt",sep=""))

# load tau0_est
tau0.est <- read.table(paste(folder,"Inference_K",K,"/tau0_est.txt",sep=""))

# load PPI_est
PPI.est <- read.table(paste(folder,"Inference_K",K,"/PPI_est.txt",sep=""))
D.est <- unlist(read.table(paste(folder,"Inference_K",K,"/IG_est.txt",sep="")))

# load X_imputed
x_imputed <- read.table(paste(folder,"Inference_K",K,"/x_imputed.txt",sep=""))


# ################################################
# # Prepare intrinsic genes for Pathway analysis #
# ################################################
# intri_gene <- gene_list[D.est==1]
# write.table(intri_gene, file = paste0("intrinsic_gene_list.txt"), row.names = FALSE,col.names = FALSE,quote = FALSE)

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
write.table(x_corrected, file = paste0(folder,"x_corrected.txt"), row.names = FALSE, col.names = FALSE)
end_time<-Sys.time()
running_time<-end_time-start_time
print(running_time)

#######
# ARI #
#######
metadata$CellType <- paste0("(",metadata$Num_Cell1,",",
                            metadata$Num_Cell2,",",
                            metadata$Num_Cell3,")")
ARI_BUSseq <- adjustedRandIndex(metadata$CellType,w_BUSseq)

########################
# Set coloring to plot #
########################

#####set cell type colorings
color_by_celltype <- c(brewer.pal(K,"Set3"))
celltype_factor <- metadata$CellType

#####set batch colorings
color_by_batch <- c("#EB4334","#FBBD06","#35AA53","#4586F3")
batch_factor <- factor(metadata$Batch)

########################################################################
# Draw t-SNE plots by batch and true cell types for the raw count data #
########################################################################
set.seed(123)
all.dists.uncorrected <- as.matrix(dist(t(log1p(y_obs))))
tsne_uncorrected_dist <- Rtsne(all.dists.uncorrected, is_distance=TRUE, perplexity = 30)

uncorrected_by_celltype<- paste0("Image/tsne_",proj,"_v",ver,"_uncorrected_by_celltype.jpeg")
dat_frame <- data.frame(Var1 = tsne_uncorrected_dist$Y[,1], 
                        Var2 = tsne_uncorrected_dist$Y[,2], 
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

uncorrected_by_batch <- paste0("Image/tsne_",proj,"_v",ver,"_uncorrected_by_batch.jpeg")
dat_frame <- data.frame(Var1 = tsne_uncorrected_dist$Y[,1], 
                        Var2 = tsne_uncorrected_dist$Y[,2], 
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

#####################
# Draw the PCA plot #
#####################
pca.uncorrected <- prcomp(t(log1p(y_obs)), rank=2)

uncorrected_PCA_by_celltype<- paste0("Image/pca_",proj,"_v",ver,"_uncorrected_by_celltype.jpeg")
dat_frame <- data.frame(Var1 = pca.uncorrected$x[,1], 
                        Var2 = pca.uncorrected$x[,2], 
                        col = celltype_factor)

jpeg(uncorrected_PCA_by_celltype,width = 800, height = 600, quality = 100)
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

uncorrected_PCA_by_batch<- paste0("Image/pca_",proj,"_v",ver,"_uncorrected_by_batch.jpeg")
dat_frame <- data.frame(Var1 = pca.uncorrected$x[,1], 
                        Var2 = pca.uncorrected$x[,2], 
                        col = batch_factor)

jpeg(uncorrected_PCA_by_batch,width = 800, height = 600, quality = 100)
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

##############################################################################
# Draw t-SNE plots by batch and true cell types for the corrected count data #
##############################################################################
set.seed(123)
all.dists.BUSseq <- as.matrix(dist(t(log1p(x_corrected[D.est==1,]))))
tsne_BUSseq_dist <- Rtsne(all.dists.BUSseq, is_distance=TRUE, perplexity = 30)

BUSseq_by_celltype<- paste0("Image/tsne_",proj,"_v",ver,"_BUSseq_by_celltype.jpeg")
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

BUSseq_by_batch <- paste0("Image/tsne_",proj,"_v",ver,"_BUSseq_by_batch.jpeg")
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

#####################
# Draw the PCA plot #
#####################
pca.BUSseq <- prcomp(t(log1p(x_corrected[D.est==1,])), rank=2)

BUSseq_PCA_by_celltype<- paste0("Image/pca_",proj,"_v",ver,"_BUSseq_by_celltype.jpeg")
dat_frame <- data.frame(Var1 = pca.BUSseq$x[,1], 
                        Var2 = pca.BUSseq$x[,2], 
                        col = celltype_factor)

jpeg(BUSseq_PCA_by_celltype,width = 800, height = 600, quality = 100)
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

BUSseq_PCA_by_batch<- paste0("Image/pca_",proj,"_v",ver,"_BUSseq_by_batch.jpeg")
dat_frame <- data.frame(Var1 = pca.BUSseq$x[,1], 
                        Var2 = pca.BUSseq$x[,2], 
                        col = batch_factor)

jpeg(BUSseq_PCA_by_batch,width = 800, height = 600, quality = 100)
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
if(ver == 2){
  jpeg(paste0("Image/legend_",proj,"_v",ver,".jpg"),width = 800, height = 600, quality = 100)
  plot(c(-8,8),c(-6,6),type="n", bty="n", axes=FALSE, xlab="", ylab="")
  legend(x=-8, y=6, legend=c("(0,0,9)", "(0,9,0)", "(9,0,0)"), pch=20, cex=2.5, pt.cex = 6, col=color_by_celltype , title = "Cell Type", bty="n")
  legend(x=0, y=6, legend=c("Mix1","Mix2","Mix3","Mix4"), pch=20, cex=2.5, pt.cex = 6, col = color_by_batch, title = "Batch", bty="n")
  dev.off()
}else{
  jpeg(paste0("Image/legend_",proj,"_v",ver,".jpg"),width = 800, height = 600, quality = 100)
  plot(c(-8,8),c(-6,6),type="n", bty="n", axes=FALSE, xlab="", ylab="")
  legend(x=-8, y=6, legend=c("(0,0,9)", "(0,9,0)","(4,0,5)","(5,0,4)", "(9,0,0)"), pch=20, cex=2.5, pt.cex = 6, col=color_by_celltype , title = "Cell Type", bty="n")
  legend(x=0, y=6, legend=c("Mix1","Mix2","Mix3","Mix4"), pch=20, cex=2.5, pt.cex = 6, col = color_by_batch, title = "Batch", bty="n")
  dev.off()
}

}
