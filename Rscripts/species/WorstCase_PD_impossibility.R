# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/species")
setwd(wd)


# Functions ---------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))

source("../../R/Rfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")


# Custom functions --------------------------------------------------------

# Plot options ------------------------------------------------------------------
save_img = FALSE
width = 12; height = 6
cex.labels <- cex.lab <- 2
cex.axis <- 2
cex.legend <- 1.5
mycol = c("darkorange","darkred","darkblue","lightblue","aquamarine")
mycol2 = c("black","lightblue")
lgd_names = c("Freq","PD","FDP","Dir-Multi")

# Options -----------------------------------------------------------------
alpha = 0.05
n = 500 
Nexp <- Nrep <- 100
M = 200
M_max = 200
Rmax = 100

seed0 <- seed <- 777
set.seed(seed0)
seeds = sample(1:999999, size = Nexp)


name = "WorstUnif"
params = c(n,alpha)
ptrue = sim_generic_species(name,M,params)
ptrue = sort(ptrue, decreasing = TRUE)
ptrue_mat = matrix(nrow = Nrep,ncol = M)
ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)
data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
Nj_mat = apply(data_mat, 2, tabulate, nbins = M) # M x Nrep
Mworst = 1/ptrue[1]

Mmax_plot = 1000
xpos = seq(1,Mmax_plot,length.out = 10); xlabs = round(xpos,0)
ymax = max(Nj_mat/n);ymin = 0
ypos = seq(ymin,ymax,length.out = 5); ylabs = round(1000*ypos,0)
ylim_plot = c(ymin,ymax)
ylab = "Emp. distr."

par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(0,(Mmax_plot+1)), ylim = ylim_plot,
     xlab = "r", ylab = ylab, 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
grid(lty = 1,lwd = 1, col = "gray90" )
for(m in 1:M){
  points(y = Nj_mat[m,]/n, x = rep(m,Nrep), pch = 16, cex = 0.5)
}
points(x = 1:M, y = ptrue, col = "red", pch = 16, type = "b", lwd = 3)

# Empirical Bayes -----------------------------------------------------------------

EB_fit = t(apply(data_mat, 2, function(data){
  n_i = tabulate(data, nbins = M)
  idx_obs = which(n_i > 0)
  Kn = length(idx_obs)
  data_obs = n_i[idx_obs]
  pmin_obs = min(ptrue[idx_obs])
  res = matrix(nrow = 1, ncol = 12)
  colnames(res) = c("sigma-PD","theta-PD",
                    "Lambda-FDP","gamma-FDP",
                    "gamma-DirMulti","Kn",
                    "Mmax","Freq",
                    "UB_PD","UB_FDP","UB_DM","pmin_obs")
  res[1,6] = Kn
  res[1,12] = pmin_obs
  #a) PD
  # Param. estimation (PYP)
  start_params <- c(alpha = 0.5, theta = 1)
  fit <- optim(par = start_params, fn = llik_pyp, 
               n = n, Kn = Kn, data_obs = data_obs, # extra parameters
               method = "L-BFGS-B",
               lower = c(0, -1), upper = c(1-1e-10, Inf)) 
  res[1,1] = fit$par[1]
  res[1,2] = fit$par[2]
  
  #b) FDP
  # Param. estimation (FDP)
  start_params <- c(gamma = 0.1, Lambda = Kn)
  fit <- optim(par = start_params, fn = llik_FD, 
               n = n, Kn = Kn, data_obs = data_obs, M_max = M_max,# extra parameters
               method = "L-BFGS-B",
               lower = c(1e-5, 1e-5), upper = c(Inf, Inf)) 
  res[1,3] = fit$par[2]
  res[1,4] = fit$par[1]  
  
  #c) Dirichlet-Multinomial --> M = M
  # Param. estimation (DirMulti)
  start_params <- c(gamma = 0.5)
  fit <- optim(par = start_params, fn = llik_DirMult,
               n = n, M = M, data = n_i, # extra parameters
               method = "L-BFGS-B",
               lower = c(1e-10), upper = c(Inf))
  res[1,5] = fit$par[1]   
  
  
  ### Calcolo Upper bounds
  if(Kn == M){
    Mmax = 0
  }else{
    idx_unobs = which(n_i == 0)
    Mmax = max(ptrue[idx_unobs])
  }
  res[1,7] = Mmax
  #b) Freq
  pain = ub_pain(n = n, Rmax = Rmax, alfa = alpha)
  pain = min(pain,1)
  res[1,8] = pain
  # Upper bound (PYP)
  ubpyp = exp(compute_log_UBMarkov( Rmax, res[1,1], res[1,2], Kn, n, alpha ))
  ubpyp = min(ubpyp,1)
  res[1,9] = ubpyp
  #d) FDP
  ubFD = exp(compute_log_UBMarkov_FD( Rmax, res[1,4], res[1,3], Kn, n, alpha, M_max ))
  ubFD = min(ubFD,1)
  res[1,10] = ubFD
  #f) Dirichlet-Multinomial --> M = M_DM
  ubDM = exp(compute_log_UB_DirMulti( Rmax, res[1,5], M, Kn, n, alpha ))
  ubDM = min(ubDM,1)
  res[1,11] = ubDM
  res
})) 
colnames(EB_fit) = c("sigma","theta","LambdaFDP","gammaFDP","gammaDM","Kn","Mmax","Freq","UB_PD","UB_FDP","UB_DM","pmin_obs")
EB_fit[c(7:11)] = EB_fit[c(7:11)] * 1e3

colMeans(EB_fit[,c(1,2)])
EB_fit[,9]

apply(EB_fit[,c(7,9,10,11)], 2, quantile, probs = c(0.025,0.5,0.975))

# MAP -----------------------------------------------------------------
var_prior = 1
MAP_fit = t(apply(data_mat, 2, function(data){
  n_i = tabulate(data, nbins = M)
  idx_obs = which(n_i > 0)
  Kn = length(idx_obs)
  data_obs = n_i[idx_obs]
  pmin_obs = min(ptrue[idx_obs])
  res = matrix(nrow = 1, ncol = 12)
  colnames(res) = c("sigma-PD","theta-PD",
                    "Lambda-FDP","gamma-FDP",
                    "gamma-DirMulti","Kn",
                    "Mmax","Freq",
                    "UB_PD","UB_FDP","UB_DM","pmin_obs")
  res[1,6] = Kn
  res[1,12] = pmin_obs
  
  ## i) PD hyperparameters
  a_sigma = 1
  b_sigma = 1
  mu_theta = Hychoice_MCMC_general("PD", n, Kn, hy_prior = c(a_sigma,b_sigma) )
  a_theta = mu_theta*mu_theta/var_prior
  b_theta = mu_theta/var_prior
  hy_pyp = c(a_sigma,b_sigma, a_theta, b_theta)
  
  ## ii) FDP hyperparameters
  mu_Lambda = Kn 
  a_Lambda = mu_Lambda*mu_Lambda/var_prior
  b_Lambda = mu_Lambda/var_prior
  mu_gamma = 1 #Hychoice_MCMC_general("FDP", n, Kn, hy_prior = c(mu_Lambda) )
  a_gamma  = mu_gamma*mu_gamma/var_prior
  b_gamma  = mu_gamma/var_prior
  hy_FDP = c(a_gamma,b_gamma, a_Lambda, b_Lambda)
  
  ## iii) DirMulti hyperparameters
  mu_gammaDM = 1
  a_gammaDM  = mu_gammaDM*mu_gammaDM/var_prior
  b_gammaDM  = mu_gammaDM/var_prior
  hy_DM = c(a_gammaDM,a_gammaDM)
  
  #a) PD
  # Param. estimation (PYP)
  start_params <- c(alpha = 0.5, theta = 1)
  fit <- optim(par = start_params, fn = lpost_pyp, 
               n = n, Kn = Kn, data_obs = data_obs, hy = hy_pyp, # extra parameters
               method = "L-BFGS-B",
               lower = c(0, -1), upper = c(1-1e-10, Inf)) 
  res[1,1] = fit$par[1]
  res[1,2] = fit$par[2]
  
  #b) FDP
  # Param. estimation (FDP)
  start_params <- c(gamma = 0.1, Lambda = Kn)
  fit <- optim(par = start_params, fn = lpost_FD, 
               n = n, Kn = Kn, data_obs = data_obs, M_max = M_max, hy = hy_FDP,  # extra parameters
               method = "L-BFGS-B",
               lower = c(1e-5, 1e-5), upper = c(Inf, Inf)) 
  res[1,3] = fit$par[2]
  res[1,4] = fit$par[1]  
  
  #c) Dirichlet-Multinomial --> M = M
  # Param. estimation (DirMulti)
  start_params <- c(gamma = 0.5)
  fit <- optim(par = start_params, fn = lpost_DirMult,
               n = n, M = M, data = n_i, hy = hy_DM, # extra parameters
               method = "L-BFGS-B",
               lower = c(1e-10), upper = c(Inf))
  res[1,5] = fit$par[1]   
  
  
  ### Calcolo Upper bounds
  if(Kn == M){
    Mmax = 0
  }else{
    idx_unobs = which(n_i == 0)
    Mmax = max(ptrue[idx_unobs])
  }
  res[1,7] = Mmax
  #b) Freq
  pain = ub_pain(n = n, Rmax = Rmax, alfa = alpha)
  pain = min(pain,1)
  res[1,8] = pain
  # Upper bound (PYP)
  ubpyp = exp(compute_log_UBMarkov( Rmax, res[1,1], res[1,2], Kn, n, alpha ))
  ubpyp = min(ubpyp,1)
  res[1,9] = ubpyp
  #d) FDP
  ubFD = exp(compute_log_UBMarkov_FD( Rmax, res[1,4], res[1,3], Kn, n, alpha, M_max ))
  ubFD = min(ubFD,1)
  res[1,10] = ubFD
  #f) Dirichlet-Multinomial --> M = M_DM
  ubDM = exp(compute_log_UB_DirMulti( Rmax, res[1,5], M, Kn, n, alpha ))
  ubDM = min(ubDM,1)
  res[1,11] = ubDM
  res
})) 
colnames(MAP_fit) = c("sigma","theta","LambdaFDP","gammaFDP","gammaDM","Kn","Mmax","Freq","UB_PD","UB_FDP","UB_DM","pmin_obs")
MAP_fit[,c(7:11)] = MAP_fit[,c(7:11)] * 1e3

colMeans(MAP_fit[,c(1,2)])

apply(MAP_fit[,c(7,9,10,11)], 2, quantile, probs = c(0.05,0.5,0.95))



# Grid of values ----------------------------------------------------------


sigma_grid = c(0,0.1,0.5,0.9)
theta_grid = seq(1e-5, 1000, length.out = 5000)
UPD = matrix(-1,nrow = length(sigma_grid), ncol = length(theta_grid))
for(i in seq_along(sigma_grid)){
  for(j in seq_along(theta_grid)){
    UPD[i,j] = exp(compute_log_UBMarkov( Rmax, sigma_grid[i], theta_grid[j], Kn=Mworst, n, alpha ))
  }
}
min(UPD)
max(UPD)

par(mfrow = c(2,2))
for(hh in 1:4){
  plot(x = theta_grid, y = UPD[hh,], type = "l", lwd = 3, ylim = c(0,1/Mworst + 1/10000 ))
  abline(h = 1/Mworst, lty = 2, col = "red", lwd = 2)
}

theta_grid[which.max(UPD[1,])]



sum(1/(72 + n + 0:10))
1/72
