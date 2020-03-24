# calculate the acceptance rate in the simulation study
rm(list=ls())
library(bigmemory)
library(xtable)

setwd("ConvergenceDiagnostic")
proj <- "hemat"
ver <- 1
K <- 6

# load the dimension information
dim <- read.table(paste0("../RawCountData/dim_",proj,"_v",ver,".txt"))
dim <- unlist(dim)
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3+1:B]

samdir <- paste0("../MCMC_sampling_K",K,"/")
setwd(samdir)

accpt_rate <- function(sampling){
  
  num_iter <- length(sampling)
  rej_ind <- sampling[2:num_iter] == sampling[1:(num_iter-1)]
  acc_rat <- 1 - sum(rej_ind)/(num_iter-1)
  return(acc_rat)
}

rate_summary <- matrix(NA, 6, 6) 

# alpha
alpha.post <- as.matrix(read.big.matrix("alpha_post.txt",sep=" ",type="double"))
alpha_rate <- apply(alpha.post, 2, accpt_rate)
rate_summary[1,] <- summary(alpha_rate) 

# beta
beta.post <- as.matrix(read.big.matrix("beta_post.txt",sep=" ",type="double"))
beta.post <- beta.post[,-(1 : G * K - K + 1)]
beta_rate <- apply(beta.post, 2, accpt_rate)
rate_summary[2,] <- summary(beta_rate) 

# nu
nu.post <- as.matrix(read.big.matrix("nu_post.txt",sep=" ",type="double"))
nu.post <- nu.post[,-(1 : G * B - B + 1)]
nu_rate <- apply(nu.post, 2, accpt_rate)
rate_summary[3,] <- summary(nu_rate) 

# delta
delta.post <- as.matrix(read.big.matrix("delta_post.txt",sep=" ",type="double"))
ref_ind <- 1
for(b in 2:B){
  ref_ind <- c(ref_ind, sum(nb[1:(b-1)])+1)
}
delta.post <- delta.post[,-ref_ind]
delta_rate <- apply(delta.post, 2, accpt_rate)
rate_summary[4,] <- summary(delta_rate) 

# phi
phi.post <- as.matrix(read.big.matrix("phi_post.txt",sep=" ",type="double"))
phi.post <- phi.post[,-(1 : G * B - B + 1)]
phi_rate <- apply(phi.post, 2, accpt_rate)
rate_summary[5,] <- summary(phi_rate)

# gamma
gamma.post <- as.matrix(read.big.matrix("gamma_post.txt",sep=" ",type="double"))
gamma_rate <- apply(gamma.post, 2, accpt_rate)
rate_summary[6,] <- summary(gamma_rate)

rownames(rate_summary) <- c("alpha","beta","nu","delta","phi","gamma")
colnames(rate_summary) <- c("Min","First Quantile","Median","Mean","Third quartile","Max")

print(rate_summary)
