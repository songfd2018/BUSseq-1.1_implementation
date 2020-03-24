# calculate the EPSR for the convergence test
rm(list=ls())
library(bigmemory) # for loading the posterior sampling
library(mclust) # for matching switched labeling
proj <- "pancreas"
ver <- 1

setwd("ConvergenceDiagnostic")

# Calculate EPSR
cal_EPSR <- function(chain1, chain2){ #EPSR for two chains
  
  num.chain <- 4
  num.iter <- nrow(chain1) / 2
  
  dim <- ncol(chain1)
  the_mean <- matrix(NA, num.chain, dim)
  the_mean[1,] <- apply(chain1[1:num.iter,],2,mean)
  the_mean[2,] <- apply(chain1[1:num.iter + num.iter,],2,mean)
  the_mean[3,] <- apply(chain2[1:num.iter,],2,mean)
  the_mean[4,] <- apply(chain2[1:num.iter + num.iter,],2,mean)
  
  the_var <- matrix(NA, num.chain, dim)
  the_var[1,] <- apply(chain1[1:num.iter,],2,var)
  the_var[2,] <- apply(chain1[1:num.iter + num.iter,],2,var)
  the_var[3,] <- apply(chain2[1:num.iter,],2,var)
  the_var[4,] <- apply(chain2[1:num.iter + num.iter,],2,var)
  
  B <- apply(the_mean, 2, var) * num.iter
  W <- apply(the_var, 2, mean) 
  EPSR <- sqrt((num.iter-1)/num.iter + B / W / num.iter)
  
  return(EPSR)
}

# set the iteration information
num.iter <- 8000
num.burnin <- 4000
m <- 4
n <- (num.iter - num.burnin) / 2

# load the dimension of count data
dim_data <- unlist(read.table(paste0("../RawCountData/dim_",proj,"_v",ver,".txt")))
N <- dim_data[1]
G <- dim_data[2]
B <- dim_data[3]
nb <- dim_data[3 + 1:B]

K <- 8

####################################
# The EPSR of logmu = alpha + beta #
####################################
message("Loading alpha and beta...")

# consider the label switching
w_est_r1 <- read.table(paste0("../Inference_K",K,"/w_est.txt"))
w_est_r2 <- read.table(paste0("Inference_K",K,"/w_est.txt"))

label_switch <- c(1,6,3,2,5,4,7,8)
w_est_r2_switched <- factor(unlist(w_est_r2 + 1), levels = label_switch)
table(w_est_r1,w_est_r2_switched)

# rep 1
variable <- "alpha"
alpha_post_r1 <- read.big.matrix(paste0("../MCMC_sampling_K",K,"/",variable,"_post.txt"),sep = " ", skip = num.burnin, type ="double")
variable <- "beta"
beta_post_r1 <- read.big.matrix(paste0("../MCMC_sampling_K",K,"/",variable,"_post.txt"),sep = " ", skip = num.burnin, type ="double")

the_post_r1 <- beta_post_r1

for(g in 1:G){
  the_post_r1[,1:K + (g-1) * K] <- the_post_r1[,1:K + (g-1) * K] + alpha_post_r1[,g]
}

# rep 2
variable <- "alpha"
alpha_post_r2 <- read.big.matrix(paste0("MCMC_sampling_K",K,"/",variable,"_post.txt"),sep = " ", skip = num.burnin, type ="double")
variable <- "beta"
beta_post_r2 <- read.big.matrix(paste0("MCMC_sampling_K",K,"/",variable,"_post.txt"),sep = " ", skip = num.burnin, type ="double")

the_post_r2 <- beta_post_r2
for(g in 1:G){
  the_post_r2[,1:K + (g-1) * K] <- the_post_r2[,1:K + (g-1) * K] + alpha_post_r2[,g]
  the_post_r2[,1:K + (g-1) * K] <- the_post_r2[,label_switch + (g-1) * K]
}
 
EPSR_mu <- cal_EPSR(the_post_r1, the_post_r2)
mean(EPSR_mu < 1.3)

##################
# The EPSR of nu #
##################
message("Loading nu...")

variable <- "nu"
# rep 1
the_post_r1 <- read.big.matrix(paste0("../MCMC_sampling_K",K,"/",variable,"_post.txt"),sep = " ", skip = num.burnin, type ="double")
num.batch <- ncol(the_post_r1) / G
the_post_r1 <- the_post_r1[,-(0:(G-1) * num.batch + 1)]

# rep 2
the_post_r2 <- read.big.matrix(paste0("MCMC_sampling_K",K,"/",variable,"_post.txt"),sep = " ", skip = num.burnin, type ="double")
the_post_r2 <- the_post_r2[,-(0:(G-1) * num.batch + 1)]

EPSR_nu <- cal_EPSR(the_post_r1, the_post_r2)
mean(EPSR_nu < 1.3)

###################
# The EPSR of phi #
###################
message("Loading phi...")

variable <- "phi"
# rep 1
the_post_r1 <- read.big.matrix(paste0("../MCMC_sampling_K",K,"/",variable,"_post.txt"),sep = " ", skip = num.burnin, type ="double")

# rep 2
the_post_r2 <- read.big.matrix(paste0("MCMC_sampling_K",K,"/",variable,"_post.txt"),sep = " ", skip = num.burnin, type ="double")

EPSR_phi <- cal_EPSR(the_post_r1, the_post_r2)
mean(EPSR_phi < 1.3)

save(EPSR_mu, EPSR_nu, EPSR_phi, file = "EPSR_factors_pancreas.RData")