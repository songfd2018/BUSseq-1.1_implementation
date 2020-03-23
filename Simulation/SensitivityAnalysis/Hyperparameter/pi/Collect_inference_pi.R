rm(list=ls())
# collect the posterior sampling results of different values of hyperparameters

library(RColorBrewer)
library(gplots)
library("devtools")
library(Rtsne)
library(ggplot2)
library(mclust)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

proj <- "simulation"
ver <- 1
par <- "pi"
hypar_values <- c(0.1,0.5,2,5,10)

# set work directory
workdir <- paste0("SensitivityAnalysis/Hyperparameter/",par,"/")
setwd(workdir)

# create the directory to store images
image_dir <- paste(workdir,"/image",sep="")
if(!dir.exists(image_dir)){
  dir.create(image_dir)
}

# create the directory to put the posterior inference
num_values<-4
k_sel <- 5
N_k <- length(k_sel)
for(k in k_sel){
  for(i in 1:num_values){
    folder_name <- paste0(par,i)
    if(!dir.exists(folder_name)){
      dir.create(folder_name)
    }
  }
}

# Load the cell type labels
metadata <- read.table(paste0("../RawCountData/metadata_",proj,"_v",ver,".txt"))
w_true <- factor(metadata$celltype)

# Load dimension
dim <- read.table(paste0("../RawCountData/dim_",proj,"_v",ver,".txt"))
dim <- unlist(dim)
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3+1:B]

######################################################
#Likelihood comparison among different initial values#
######################################################
BIC.record <- rep(NA,num_values)
num_intri <- rep(NA,num_values)
tau0sq_record <- rep(NA,num_values)
ARI <- rep(NA,num_values)

#dropout_rate <- array(NA, dim=c(N_k,rep,B)) # dropout rate per batch
names(BIC.record) <- hypar_values[-3]
names(num_intri) <- hypar_values[-3]
names(tau0sq_record) <- hypar_values[-3]
names(tau0sq_record) <- hypar_values[-3]
names(ARI) <- hypar_values[-3]


for(i in 1:num_values){
    post_dir<-paste0(workdir,par,i,"/Inference_K",k_sel,"/")
    if(dir.exists(post_dir)){
      # load BIC
      BIC<- unlist(read.table(paste0(post_dir,"BIC.txt")))
      BIC.record[i] <- BIC[1]
        
      # load intrinsic genes
      IG_est<- unlist(read.table(paste0(post_dir,"IG_est.txt")))
      num_intri[i] <- sum(IG_est)

      # load tau0_est
      tau0sq_record[i]<- unlist(read.table(paste0(post_dir,"tau0_est.txt")))

      # load cell type indicators and calculate ARI
      w_est <- unlist(read.table(paste0(post_dir,"w_est.txt")))
      ARI[i]<- adjustedRandIndex(w_true,w_est)
        
      message(paste0("Finish loading the results of ",i,"-th value of ",par,"!"))
    }else{
      message(paste("The results of ",i,"-th value of ",par," doesn't exist!/n",sep=""))
    }
}

# Export to Excel by clipboard
Output <- matrix(NA, 3, num_values)
rownames(Output) <- c("BIC","ARI","No. intrinsic genes")
colnames(Output) <- hypar_values[-3]
Output[1,] <- BIC.record
Output[2,] <- ARI
Output[3,] <- num_intri

print(Output)



# draw scatter plot of ARI
ARI_values <- c(ARI[1:2],1,ARI[3:4])

jpeg(paste0(image_dir,"/Sensitivity_analysis_",proj,"_v",ver,"_",par,".jpg"),width = 540, height = 720)
par(mar=c(5.1,5.1,2.1,2.1))
plot(1:5,ARI_values,
     xlab="Hyperparameter values", xaxt = "none",
     ylab = "ARI",ylim = c(0,1), 
     type="l", cex.axis = 2, cex.lab = 3)
points(1:5,ARI_values,cex=3,pch = 19)
axis(1, at = 1:5, hypar_values, cex.axis = 2)
dev.off()