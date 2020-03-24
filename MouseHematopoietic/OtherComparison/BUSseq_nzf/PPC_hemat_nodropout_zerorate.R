######################################
# Conduct posterior prediction check #
######################################
set.seed(6543)
rm(list=ls())
proj <- "hemat"
K <- 6
ver <- 1

setwd("OtherComparison/BUSseq_nzf")

# load dimension information
dim <- unlist(read.table(paste0("../../RawCountData/dim_",proj,"_v",ver,".txt")))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# open all files
# pi 
file_pi=paste0("MCMC_sampling_K",K,"/pi_post.txt")
con_pi=file(file_pi,open="r")

# alpha
file_alpha=paste0("MCMC_sampling_K",K,"/alpha_post.txt")
con_alpha=file(file_alpha,open="r")

# beta
file_beta=paste0("MCMC_sampling_K",K,"/beta_post.txt")
con_beta=file(file_beta,open="r")

# nu
file_nu=paste0("MCMC_sampling_K",K,"/nu_post.txt")
con_nu=file(file_nu,open="r")

# delta
file_delta=paste0("MCMC_sampling_K",K,"/delta_post.txt")
con_delta=file(file_delta,open="r")

# phi
file_phi=paste0("MCMC_sampling_K",K,"/phi_post.txt")
con_phi=file(file_phi,open="r")

# simulate data for each posterior sampling after the burn-in iterations
tot_iter <- 8000
burnin <- 4000
# skip burn-in
start_skip_burnin <- Sys.time()
for(i in 1:burnin){
  pi.iter <- readLines(con_pi, n = 1)
  # w.iter <- readLines(con_w, n = 1)
  alpha.iter <- readLines(con_alpha, n = 1)
  beta.iter <- readLines(con_beta, n = 1)
  nu.iter <- readLines(con_nu, n = 1)
  delta.iter <- readLines(con_delta, n = 1)
  phi.iter <- readLines(con_phi, n = 1)
}
end_skip_burnin <- Sys.time()
print(paste0("It takes ",difftime(end_skip_burnin, start_skip_burnin, units = "mins")," mins to skip burnin iterations."))

# record the zero rate of each replication
num.output <- 100
zero_rate <- matrix(NA, num.output, B + 1)

y.temp <- matrix(NA,G,N)
z.temp <- matrix(NA,G,N)
x.temp <- matrix(NA,G,N)
w.temp <- matrix(NA,N)

start_simulate <- Sys.time()
iter <- 1
for(i in 1:(tot_iter - burnin)){
# for(i in 1:100){ 
  # read data from file
  pi.iter <- readLines(con_pi, n = 1)
  alpha.iter <- readLines(con_alpha, n = 1)
  beta.iter <- readLines(con_beta, n = 1)
  nu.iter <- readLines(con_nu, n = 1)
  delta.iter <- readLines(con_delta, n = 1)
  phi.iter <- readLines(con_phi, n = 1)
  
  # change character lines to vectors or matrices
  pi.temp <- matrix(as.numeric(unlist(strsplit(pi.iter, split = " "))), K, B)
  alpha.temp <- as.numeric(unlist(strsplit(alpha.iter, split = " ")))
  beta.temp <- matrix(as.numeric(unlist(strsplit(beta.iter, split = " "))), G, K)
  nu.temp <- matrix(as.numeric(unlist(strsplit(nu.iter, split = " "))), G, B)
  delta.temp <- as.numeric(unlist(strsplit(delta.iter, split = " ")))
  phi.temp <- matrix(as.numeric(unlist(strsplit(phi.iter, split = " "))), G, B)
  
  # simulate the observation replications
  ind_cell <- 1
  for(b in 1:B){
    for(i in 1:nb[b]){
      # w_bi ~ multinomial(pi_b)
      w.temp[ind_cell] <- which(rmultinom(1, size = 1, prob = pi.temp[,b])==1)
      
      # mu_big = exp( alpha_g + beta_gk + nu_bg + delta_bi )
      mu.temp <- exp(alpha.temp + beta.temp[,w.temp[ind_cell]] + nu.temp[,b] + delta.temp[ind_cell])
      
      # y_big ~ NB(mu_big, phi_bg)
      y.temp[,ind_cell] <- rnbinom(G, size = phi.temp[,b], mu = mu.temp)
      
      ind_cell <- ind_cell + 1      
    }
  }
  
  # Calculate the zero rate in each replication
  zero_rate[iter, 1] <- sum(y.temp == 0) / N / G
  
  ind_first <- 0
  for(b in 1:B){
    
    ind_batch <- ind_first + 1:nb[b]
    zero_rate[iter, b + 1] <- sum(y.temp[,ind_batch] == 0)/ nb[b] / G
    ind_first <- ind_first + nb[b]
    
  }
  
  if(iter == num.output){
    end_simulate <- Sys.time()
    print(paste0("It takes ",difftime(end_simulate, start_simulate, units = "mins")," mins to simulate ",iter," replications."))
    write.table(zero_rate, file = paste0("ZeroRate_",proj,"_nzf.txt"), append = TRUE, row.names = FALSE, col.names = FALSE)
    iter <- 0
  }
  iter <- iter + 1
}

close(con_pi)
close(con_alpha)
close(con_beta)
close(con_nu)
close(con_delta)
close(con_phi)