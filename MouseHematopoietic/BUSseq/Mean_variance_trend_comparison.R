# Compare the mean-variance trend of real data with that of simulated data by BUSseq model
rm(list=ls())
library(ggplot2)

setwd("BUSseq")
proj <- "Hemat"
b <- 2

load("../BUSseq_workspace.RData")

set.seed(123)
Num_sel <- table(w_BUSseq[1:nb[2]+nb[1]])
# Consider the second batch
cell_sel <- list()
for(k in 1:K){
  
  cell_sel[[k]] <- which(w_BUSseq[1:nb[2]+nb[1]] == k - 1)
}

# true expression levels 
true_count <- list()

# cell library size
delta_sel <- list()

for(k in 1:K){
  true_count[[k]] <- y_obs[,cell_sel[[k]]]
  delta_sel[[k]] <- delta.est[cell_sel[[k]]]
}


# simulate the expression levels of cells
simulate_count <- list()


for(k in 1:K){
  y_temp <- matrix(NA, G, Num_sel[k])
  for(i in 1:Num_sel[k]){
    for(g in 1:G){
      # underlying expression levels
      mean_level <- exp(alpha.est[g] + beta.est[g,k] + nu.est[g, b] + delta_sel[[k]][i])
      x <- rnbinom(1,phi.est[g,b], mu= mean_level)
      
      # dropout events
      logit_pi <- gamma.est[b,1] + gamma.est[b,2] * x
      drop_ind <- rbinom(1,1,prob = 1/(1+exp(-logit_pi)))
      if(drop_ind==0){
        y_temp[g,i] <- x
      }else{
        y_temp[g,i] <- 0
      }
    }
  }
  simulate_count[[k]] <- y_temp
}

# mean-variance levels
if(!dir.exists("Image")){
  dir.create("Image")
}

log_true_count_collect <- matrix(NA, G, sum(Num_sel))
log_simulate_count_collect <- matrix(NA, G, sum(Num_sel))
temp_index <- 0
pl <- list()

for(k in 1:K){
  #Count Per Million normalization
  log_true_count_normalized <- true_count[[k]]
  log_simulate_count_normalized <- simulate_count[[k]]
  for(i in 1:Num_sel[k]){
    temp_sum <- sum(log_true_count_normalized[,i])
    log_true_count_normalized[,i] <- log1p(log_true_count_normalized[,i]/ temp_sum * 10^6)
    temp_sum <- sum(log_simulate_count_normalized[,i])
    log_simulate_count_normalized[,i] <- log1p(log_simulate_count_normalized[,i]/ temp_sum * 10^6)
  }
  
  # collect the normalized count data 
  log_true_count_collect[, temp_index + 1:Num_sel[k]] <- as.matrix(log_true_count_normalized)
  log_simulate_count_collect[, temp_index + 1:Num_sel[k]] <- log_simulate_count_normalized
  
  # calculate the gene-specific mean and variance
  mean_true <- apply(log_true_count_normalized,1,mean)
  var_true <- apply(log_true_count_normalized,1,var)
  
  mean_simulate <- apply(log_simulate_count_normalized,1,mean)
  var_simulate <- apply(log_simulate_count_normalized,1,var)
  
  # mean-variance trend, first plot simulated data, then plot real data
  meanvar_frame <- data.frame(Mean = c(mean_simulate, mean_true),
                              Var = c(var_simulate, var_true),
                              Category = rep(c("Simulated","Real"),each = G),
                              Cluster = paste0("Clster_",rep(k, 2 * G)))
  
  # reorder nodes for visualization
  reorder <- sample(nrow(meanvar_frame))
  meanvar_frame <- meanvar_frame[reorder,]
  
  pl[[k]] <- ggplot(meanvar_frame, aes(x= Mean, y= Var, colour= Category)) +
    geom_point(size=4, alpha = 0.3) + theme_classic() +
    scale_colour_manual(values = c("#F8766D","#00BFC4")) +
    xlab("Mean") + ylab("Variance") +
    ggtitle(paste0("Cluster ",k)) +
    theme(plot.title = element_text(hjust = 0.5, face="bold", size=48),
          axis.text.x = element_text(face = "bold", #color = "#993333", 
                                     size = 36), #angle = 45),
          axis.text.y = element_text(face = "bold", #color = "blue", 
                                     size = 36),#, angle = 45))
          axis.title=element_text(size=44,face="bold"),
          legend.text=element_text(size=28,face="bold"),
          legend.title = element_text(size = 36,face="bold"),
          panel.background = element_rect(colour = "black",size = 2))
  
  temp_index <- temp_index + Num_sel[k]
}

# Cluster 1
k <- 1
jpeg(paste0("Image/Mean_variance_trend_comparison_",proj,"_cluster_",k,".jpg"),width = 800, height = 600, quality = 100)
pl[[k]]
dev.off()


# Cluster 2
k <- 2
jpeg(paste0("Image/Mean_variance_trend_comparison_",proj,"_cluster_",k,".jpg"),width = 800, height = 600, quality = 100)
pl[[k]]
dev.off()


# Cluster 3
k <- 3
jpeg(paste0("Image/Mean_variance_trend_comparison_",proj,"_cluster_",k,".jpg"),width = 800, height = 600, quality = 100)
pl[[k]]
dev.off()

# Cluster 4
k <- 4
jpeg(paste0("Image/Mean_variance_trend_comparison_",proj,"_cluster_",k,".jpg"),width = 800, height = 600, quality = 100)
pl[[k]]
dev.off()


# Cluster 5
k <- 5
jpeg(paste0("Image/Mean_variance_trend_comparison_",proj,"_cluster_",k,".jpg"),width = 800, height = 600, quality = 100)
pl[[k]]
dev.off()


# Cluster 6
k <- 6
jpeg(paste0("Image/Mean_variance_trend_comparison_",proj,"_cluster_",k,".jpg"),width = 800, height = 600, quality = 100)
pl[[k]]
dev.off()



# mean-variance trend when combining all clusters
mean_true <- apply(log_true_count_collect,1,mean)
var_true <- apply(log_true_count_collect,1,var)

mean_simulate <- apply(log_simulate_count_collect,1,mean)
var_simulate <- apply(log_simulate_count_collect,1,var)

meanvar_frame <- data.frame(Mean = c(mean_simulate, mean_true),
                            Var = c(var_simulate,var_true),
                            Category = rep(c("Simulated","Real"),each = G),
                            Cluster = paste0("Combined"))

# reorder nodes for visualization
reorder <- sample(nrow(meanvar_frame))
meanvar_frame <- meanvar_frame[reorder,]

jpeg(paste0("Image/Mean_variance_trend_comparison_",proj,"_combined.jpg"),width = 800, height = 600, quality = 100)
ggplot(meanvar_frame, aes(x= Mean, y= Var, colour= Category)) +
  geom_point(size=4, alpha = 0.1) + theme_classic() +
  scale_colour_manual(values = c("#F8766D","#00BFC4")) +
  xlab("Mean") + ylab("Variance") + ggtitle("Hematopoietic") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=48),
        axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"),
        axis.title=element_text(size=44,face="bold"),
        panel.background = element_rect(colour = "black",size = 2))
dev.off()