rm(list=ls())
# heatmap.3 referring https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

library(ggplot2)
library(reshape2)

setwd("SensitivityAnalysis/ModelMisspecification")
ver <- 2

set.seed(12357)
################################
# Set the Synthetic Parameters #
################################
#The number of batches
B<-4

#The number of cells per batch
nb<-c(300,300,200,200)

#The total number of cells
N<-sum(nb)

#The number of genes
G<-3000

#The number of cell types
K<-5

#The first column of gamma.syn denotes the intercept of 
#the logistic regression for dropout events
#The second column of gamma.syn denotes the odds ratios 
#of the logistic regression for dropout events
gamma.syn<-matrix(0,B,2)
gamma.syn[1,]<-c(-0.5,-0.2)
gamma.syn[2,]<-c(-0.5,-0.2)
gamma.syn[3,]<-c(-0.5,-0.2)
gamma.syn[4,]<-c(-0.5,-0.2)

#the log-scale baseline expression levels
alpha.syn<-rep(NA,G)
alpha.syn[1:(G/5)]<-rep(1.5,G/5)
alpha.syn[(G/5+1):(G/5*2)]<-rep(1,G/5)
alpha.syn[(G/5*2+1):(G/5*3)]<-rep(0.5,G/5)
alpha.syn[(G/5*3+1):(G/5*4)]<-rep(0.25,G/5)
alpha.syn[(G/5*4+1):G]<-rep(0,G/5)

alpha.syn[1:50] <- 3.5
alpha.syn[G/5 + 1:50] <- 3
alpha.syn[G/5 * 2 + 1:25] <- 2.5
alpha.syn[G/5 * 3 + 1:25] <- 2.25

#the cell-type effects 
beta.syn<-matrix(0,G,K)

#the first cell type is regarded as the reference cell type 
#without cell-type effects 
beta.syn[,1] <- 0

#the cell-type effects of the second cell type
beta.syn[1:25,2] <- -2
beta.syn[26:50, 2] <- -1.5
beta.syn[51:75, 2] <- 2
beta.syn[G/5 + 1:25,2] <- -2
beta.syn[G/5 + 26:50,2] <- -1.5
beta.syn[G/5 + 51:75,2] <- 2 
beta.syn[G/5 * 2 + 1:25,2] <- -2
beta.syn[G/5 * 2 + 26:50,2] <- 2
beta.syn[G/5 * 3 + 1:25,2] <- -2
beta.syn[G/5 * 3 + 26:50,2] <- 2

#the cell-type effects of the third cell type
beta.syn[1:50,3] <- -2
beta.syn[51:75, 3] <- 2
beta.syn[G/5 + 1:50,3] <- -2
beta.syn[G/5 + 51:75,3] <- 2 
beta.syn[G/5 * 2 + 1:25,3] <- -2
beta.syn[G/5 * 2 + 46:75,3] <- 2 
beta.syn[G/5 * 3 + 1:25,3] <- -2
beta.syn[G/5 * 3 + 46:75,3] <- 2 

#the cell-type effects of the forth cell type
beta.syn[1:50,4] <- -2
beta.syn[76:95, 4] <- 2
beta.syn[96:100, 4] <- 1
beta.syn[G/5 + 1:50,4] <- -2
beta.syn[G/5 + 76:100,4] <- 2 
beta.syn[G/5 * 2 + 1:25,4] <- -2
beta.syn[G/5 * 2 + 76:100,4] <- 2  
beta.syn[G/5 * 3 + 1:25,4] <- -2
beta.syn[G/5 * 3 + 76:100,4] <- 2

#the cell-type effects of the fifth cell type
beta.syn[1:10,5] <- -1
beta.syn[11:50,5] <- -2
beta.syn[101:125, 5] <- 2 
beta.syn[G/5 + 1:15,5] <- -1
beta.syn[G/5 + 16:50,5] <- -2
beta.syn[G/5 + 101:125,5] <- 2
beta.syn[G/5 * 2 + 1:25,5] <- -2
beta.syn[G/5 * 2 + 86:125,5] <- 2
beta.syn[G/5 * 3 + 1:25,5] <- -2
beta.syn[G/5 * 3 + 86:125,5] <- 2 

# Check the values of alpha + beta
plot(1:G, alpha.syn,col=1,ylim = c(-1,5))
points(1:G, alpha.syn + beta.syn[, 2],col=2)
points(1:G, alpha.syn + beta.syn[, 3],col=3)
points(1:G, alpha.syn + beta.syn[, 4],col=4)
points(1:G, alpha.syn + beta.syn[, 5],col=5)

# generate the highly expressed indicators
HE_ind <- matrix(0, G, K)
logmu.syn <- alpha.syn + beta.syn
mean_level <- apply(logmu.syn, 1, mean)
for(k in 1:K){
  HE_ind[,k] <- logmu.syn[,k] > mean_level
}

#the batch effects
nu.syn<-matrix(NA,G,B)

#the first batch is taken as the reference batch 
#without batch effects
nu.syn[,1] <- 0

#the batch effect of the second batch
nu.syn[,2] <- rep(c(4,3,2,1,0),each = G/5)

#the batch effect of the third batch
nu.syn[,3] <- rep(c(2,1,0,1,2),each = G/5)

#the batch effect of the forth batch
nu.syn[,4] <- rep(c(-1,-2,-3,-2,-1),each = G/5)

#the cell-specific size factors
delta.syn <- list()
for(b in 1:B){
  delta.syn[[b]] <- rep(NA, nb[b])
}


#the first cell in each batch is regarded as the reference cell 
#with the cell-specific size factors being 0
delta.syn[[1]][1:100] <- 0
delta.syn[[1]][101:200] <- 1
delta.syn[[1]][201:300] <- -1

#the second batch
delta.syn[[2]][1:50] <- 0
delta.syn[[2]][51:70] <- -2
delta.syn[[2]][71:100] <- -1
delta.syn[[2]][101:300] <- -3

#the third batch
delta.syn[[3]][1:30] <- 0
delta.syn[[3]][31:90] <- -2
delta.syn[[3]][91:150] <- -1
delta.syn[[3]][151:200] <- -3


#the forth batch
delta.syn[[4]][1:70] <- 0
delta.syn[[4]][71:150] <- 1
delta.syn[[4]][151:200] <- 2

#the batch-specific and gene-specific overdispersion parameters
phi.syn<-matrix(5,B,G)#mean 2 var 0.5
phi.syn[,(G*0.4 + 1):(G * 0.8)]<-3
phi.syn[,(G*0.8 + 1):G]<-1

#the cell-type proportions in each batch
pi.syn <- matrix(NA,K,B)

#the first batch
pi.syn[,1]<-c(0.4,0.3,0.2,0.1,0)

#the second batch
pi.syn[,2]<-c(0,0.2,0.3,0.3,0.2)

#the third batch
pi.syn[,3]<-c(0.24,0,0.2,0.26,0.3)

#the forth batch
pi.syn[,4]<-c(0,0.3,0.4,0.3,0)


##############################################
# Simulate Latent Varibles and Observed data #
##############################################
#the cell-type indicators of each cell
w <- list()

#the first batch
w[[1]] <- rep(1:K, nb[1] * pi.syn[,1])

#the second batch
w[[2]] <- rep(1:K, nb[2] * pi.syn[,2])

#the third batch
w[[3]] <- rep(1:K, nb[3] * pi.syn[,3])

#the forth batch
w[[4]] <- rep(1:K, nb[4] * pi.syn[,4])

#the indicators for dropout events
z<-list()

#the underlying true expression levels
x<-list()

#the observed expression levels
y<-list()

#the logarithm of mean expreesion level of each gene in each cell
log.mu<-list()

for(b in 1:B){
  z[[b]] <- matrix(NA, G, nb[b])
  x[[b]] <- matrix(NA, G, nb[b])
  y[[b]] <- matrix(NA, G, nb[b])
  log.mu[[b]] <- matrix(NA, G, nb[b])
}

#generate the latent variable and observed data
for(b in 1:B){
  for(i in 1:nb[b]){
    log.mu[[b]][,i] <- alpha.syn + beta.syn[,w[[b]][i]] 
    log.mu[[b]][,i] <- log.mu[[b]][,i] + nu.syn[,b]
    log.mu[[b]][,i] <- log.mu[[b]][,i] + delta.syn[[b]][i]
    
    for(j in 1:G){
      if(HE_ind[j,w[[b]][i]] == 1){
        x[[b]][j,i]<-rnbinom(1,phi.syn[b,j] * 5,
                             mu=exp(log.mu[[b]][j,i]))
      }else{
        x[[b]][j,i]<-rnbinom(1,phi.syn[b,j],
                             mu=exp(log.mu[[b]][j,i]))
      }
      
      logit_pi <- gamma.syn[b,1] + gamma.syn[b,2] * x[[b]][j,i]
      
      z[[b]][j,i]<-rbinom(1,1,prob = 1/(1+exp(-logit_pi)))
      if(z[[b]][j,i]==0){
        y[[b]][j,i]<- x[[b]][j,i]
      }else{
        y[[b]][j,i]<- 0
      }
    }
  }
}

# Dropout rate
sum(unlist(z)==1)/N/G
for(b in 1:B){
  print(sum(unlist(z[[b]])==1)/nb[b]/G)
}


# Zero rate
sum(unlist(y)==0)/N/G
for(b in 1:B){
  print(sum(unlist(y[[b]])==0)/nb[b]/G)
}

if(!dir.exists("RawCountData")){
  dir.create("RawCountData")
}


save.image(paste0("simulation_countdata_v",ver,".RData"))

# Write out the read count matrix, dimsion information and metadata
readcount <- do.call(cbind,y)
write.table(readcount, file = paste0("RawCountData/count_data_simulation_v",ver,".txt"),row.names = FALSE,col.names = FALSE)

dim <- c(N,G,B,nb)
write.table(dim, file = paste0("RawCountData/dim_simulation_v",ver,".txt"),row.names = FALSE, col.names = FALSE)

metadata <- data.frame( batch = paste0("Batch_",rep(1:B,nb)),
                        celltype = paste0("Type_",unlist(w)))
write.table(metadata, file = paste0("RawCountData/metadata_simulation_v",ver,".txt"))



######################################################################
# Draw the Heatmap for the True Parameter Values and True Count Data #
######################################################################
if(!dir.exists("Image")){
  dir.create("Image")
}

# batch color bar
color_by_batch<-c("#EB4334","#FBBD06","#35AA53","#4586F3")
batch.cols<-color_by_batch[rep(1:B,nb)]

# cell type color bar
color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978","#8097D3","#9C7ACE")

############################
# Draw Heatmap in Figure 3 #
############################
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
break_mu <- seq(0, 3.5, length.out = 12)
break_nu <- seq(-4.2, 4.2, length.out = 12)
break_logcount <- seq(0,8.5,length.out = 12)

######################
# Plot the color key #
scale01 <- function(x, low = min(x), high = max(x)) {
  x <- (x - low)/(high - low)
  x
}
# color key for figure 3 (a,f)
z <- seq(min(break_mu), max(break_mu), length = length(color_key_mu))
jpeg("Image/heatmap_simulation_col_key_mu.jpg",width = 360, height = 270, quality =100)
image(z = matrix(z, ncol = 1), col = color_key_mu, breaks = break_mu, xaxt = "n", yaxt = "n")
lv <- pretty(break_mu)
xv <- scale01(as.numeric(lv), min(break_mu), max(break_mu))
axis(1, at = xv, labels = lv, cex.axis = 3)
mtext(side = 1, "Value", line = 3.5, cex = 4)
title(main = "Color Key", cex.main = 4)
dev.off()


# color key for figure 3 (b,g)
z <- seq(min(break_nu), max(break_nu), length = length(color_key_nu))
jpeg("Image/heatmap_simulation_col_key_nu.jpg",width = 360, height = 270, quality =100)
image(z = matrix(z, ncol = 1), col = color_key_nu, breaks = break_nu, xaxt = "n", yaxt = "n")
lv <- pretty(break_nu)
xv <- scale01(as.numeric(lv), min(break_nu), max(break_nu))
axis(1, at = xv, labels = lv, cex.axis = 3)
mtext(side = 1, "Value", line = 3.5, cex = 4)
title(main = "Color Key", cex.main = 4)
dev.off()

# color key for figure 3 (c-d,h-i)
z <- seq(min(break_logcount), max(break_logcount), length = 10)
jpeg("Image/heatmap_simulation_col_key_read_count.jpg",width = 360, height = 270, quality =100)
image(z = matrix(z, ncol = 1), col = colorsChoice_count, breaks = break_logcount, xaxt = "n", yaxt = "n")
lv <- pretty(break_logcount)
xv <- scale01(as.numeric(lv), min(break_logcount), max(break_logcount))
axis(1, at = xv, labels = lv, cex.axis = 3)
mtext(side = 1, "Value", line = 3.5, cex = 4)
title(main = "Color Key", cex.main = 4)
dev.off()

############################################################
# the heatmap of log-scale mean expression levels (Fig 3a) #
log_mu_celltype <- alpha.syn + beta.syn
jpeg("Image/heatmap_simulation_log_mean_expression_levels.jpg",width = 1080, height = 1440, quality =100)
heatmap.3(log_mu_celltype,
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = color_by_celltype,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key 
          lhei=c(0.4,0.4,3.6),
          col=color_key_mu,breaks = break_mu, key = FALSE)
dev.off()

##################################################
# the heatmap of location batch effects (Fig 3b) #

jpeg("Image/heatmap_simulation_batch_effects.jpeg",width = 800, height = 1000, quality =100)
heatmap.3(nu.syn,
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = color_by_batch,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key 
          lhei=c(0.6,0.4,3.6),
          col=color_key_nu,breaks = break_nu, key = FALSE)
dev.off()

#####################################################################
# the heatmap of log-scale underlying true read count data (Fig 3c) #
log_x <- NULL
for(b in 1:B){
  log_x <- cbind(log_x, log1p(x[[b]]))
}

#generate the double color bar
#batch color bar for count data matrix
color_by_batch_count <- rep(color_by_batch,nb)

#cell type color bar for count data matrix
color_by_celltype_count <- NULL
for(b in 1:B){
  color_by_celltype_count <- c(color_by_celltype_count, color_by_celltype[w[[b]]])
}

#upper color bar for batch, lower color bar for cell type
col_annotation<-cbind(color_by_celltype_count,color_by_batch_count)
colnames(col_annotation) <- NULL

jpeg("Image/heatmap_simulation_log_underlying_true_count.jpeg",width = 1080, height = 1440, quality =100)
heatmap.3(log_x,
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = col_annotation,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
          lhei=c(0.6,0.4,3.6),
          col=colorsChoice_count,breaks = break_logcount, key = FALSE)
dev.off()

##############################################################
# the heatmap of log-scale observed read count data (Fig 3d) #
log_y <- NULL
for(b in 1:B){
  log_y <- cbind(log_y, log1p(y[[b]]))
}

#share the same double color bar and splitting points 
#as that of the underlying true count data

jpeg("Image/heatmap_simulation_log_observed_count.jpeg",width = 1080, height = 1440, quality =100)
heatmap.3(log_y,
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = col_annotation,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
          lhei=c(0.5,0.4,3.6),
          col=colorsChoice_count,breaks = break_logcount, key = FALSE)
dev.off()


##################################################
# Save the Cell-specific Size Factors for Fig 3e #
##################################################
if(!dir.exists("True_para")){
  dir.create("True_para")
}

#write delta.syn as a column vector
delta.file <- "./True_para/delta_syn.txt"
file.create(delta.file)

for(b in 1:B){
  write.table(delta.syn[[b]],delta.file,append = T, col.names = F, row.names = F)
}

#write w as a column vector
w.file <- "./True_para/w_syn.txt"
file.create(w.file)


###########################
# jitter plot Figure 4(a) #
###########################
beta_nu_syn <- cbind(beta.syn,nu.syn)
colnames(beta_nu_syn) <- c(paste0("CellType",1:K),paste0("Batch",1:B))

# exclude the first cell type and first batch
beta_nu_syn <- beta_nu_syn[,-c(1,K+1)]

beta_nu_syn_frame <- melt(beta_nu_syn)
colnames(beta_nu_syn_frame) <- c("Gene","Var","Value")

jpeg("Image/simulation_jitter_beta_nu_comparison.jpg",width = 900, height = 600, quality =100)
p <- ggplot(beta_nu_syn_frame, aes(x = Var, y = Value, col = Var)) + 
  geom_jitter(shape = 20, size = 0.5) + theme_bw() +
  ylim(-5,5) + 
  geom_vline(xintercept = 4.5, linetype="dashed", size = 1) +
  scale_color_manual(values = c(color_by_celltype[-1],color_by_batch[-1])) +
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
p
dev.off()

for(b in 1:B){
  write.table(w[[b]],w.file,append = T, col.names = F, row.names = F)
}