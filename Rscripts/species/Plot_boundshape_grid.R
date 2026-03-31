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



# Plot options ------------------------------------------------------------------
save_img = FALSE
width = 12; height = 6
cex.labels <- cex.lab <- 1.5
cex.axis <- 1.5
cex.legend <- 1.5
mycol = c("darkorange","darkred","darkblue","lightblue")
mycol_var = hcl.colors(n = 4, palette = "Greens", rev = TRUE)
mycol2 = c("black","lightblue")
lgd_names = c("Oracle","Freq","PD","FDP","Dir-Multi")


# Options -----------------------------------------------------------------
n = 500
Rmax = 100
Nrep = 200
alpha <- alfa <- 0.05
seed = 42
M_max = 200

M = 200
Kn = 100
## DM univariate  -------------------------------------------------------
gamma_grid = c(seq(0.001,1,length.out = 100),
               seq(1,150,length.out = 200))
U_DM_grid = rep(NA, length(gamma_grid))
for(hh in seq_along(gamma_grid)) {
  gamma = gamma_grid[hh]
  ubFDP = exp(compute_log_UB_DirMulti( Rmax, gamma, M, Kn, n, alpha ))
  ubFDP = min(ubFDP,1)
  U_DM_grid[hh] = ubFDP
}
U_DM_grid = U_DM_grid*1000

if(save_img)
  pdf("img/UB_shape_DM.pdf", width = width, height = height)
par( mfrow = c(1,1), mar = c(2.5,3,1,1), mgp=c(1.5,0.55,0), bty = "l", las = 1, cex = 2 )
plot(gamma_grid, U_DM_grid, type="l", lwd=3,
     xlim=c(0,150), ylim=c(2,15),
     xlab="",  
     ylab=expression(U^{Dir}))
mtext(expression(gamma), side=1, line=1.25, cex = 2)  # control xlab distance here
if(save_img)
  dev.off()

## FDP - bivariate -------------------------------------------------------
Kn = 200
Lambda_grid = c(50,150,200,250,300,500)
gamma_grid = c(seq(0.001,1,length.out = 100),
               seq(1,150,length.out = 200))
res = vector("list", length = length(Lambda_grid))
counter = 1
for(gg in seq_along(Lambda_grid)){
  Lambda = Lambda_grid[gg]
  for(hh in seq_along(gamma_grid)) {
    gamma = gamma_grid[hh]
    ubFDP = exp(compute_log_UBMarkov_FD( Rmax, gamma, Lambda, Kn, n, alpha, M_max ))
    ubFDP = min(ubFDP,1)
    res[[gg]] = c(res[[gg]],ubFDP)
  }
}
res = lapply(res, function(x) x*1000)

if(save_img)
  pdf("img/UB_shape_FDP.pdf", width = width*3, height = height)
par( mfrow = c(1,3), mar = c(2.5,3,1,1), mgp=c(1.5,0.55,0), bty = "l", las = 1, cex = 2 )
for(gg in c(1,3,6)){
  plot(x = gamma_grid, y = res[[gg]], type = "l", lwd = 3 ,
       ylab = bquote(U^{FDP} ~ "- (" * Lambda == .(Lambda_grid[gg]) * ")"), 
       ylim = c(2,12.5), 
       xlab = "")
  mtext(expression(gamma), side=1, line=1.25, cex = 2)  # control xlab distance here
}
if(save_img)
  dev.off()

## PD - bivariate -------------------------------------------------------
Kn = 70 # -> la forma non cambia ma il massimo si, in funzione di Kn
sigma_grid = c(0,0.25,0.5,0.7,0.9,0.99)
theta_grid = c(seq(1e-10,1,length.out = 1000),
               seq(1,1000,length.out = 1000))
res = vector("list", length = length(sigma_grid))
counter = 1
for(gg in seq_along(sigma_grid)){
  sigma = sigma_grid[gg]
  for(hh in seq_along(theta_grid)) {
    theta = theta_grid[hh]
    ubPD = exp(compute_log_UBMarkov( Rmax, sigma, theta, Kn, n, alpha ))
    ubPD = min(ubPD,1)
    res[[gg]] = c(res[[gg]],ubPD)
  }
}
res = lapply(res, function(x) x*1000)

if(save_img)
  pdf("img/UB_shape_PD.pdf", width = 3*width, height = height)
par( mfrow = c(1,3), mar = c(2.5,3,1,1), mgp=c(1.5,0.55,0), bty = "l", las = 1, cex = 2 )
for(gg in c(1,3,5)){
  plot(x = theta_grid, y = res[[gg]], type = "l", lwd = 3 ,
       ylab = bquote(U^{PD} ~ "- (" * sigma == .(sigma_grid[gg]) * ")"), 
       ylim = c(2,12.5), 
       xlab = "")
  mtext(expression(theta), side=1, line=1.25, cex = 2)  # control xlab distance here
}
if(save_img)
  dev.off()


# Likelihood plot - DirMulti ---------------------------------------------------------

M = 200
Nrep = 50
n = 500
name = "NegBin"
params = c(1,0.003)
ptrue = sim_generic_species(name,M,params)
ptrue = sort(ptrue, decreasing = TRUE)
ptrue_mat = matrix(nrow = Nrep,ncol = M)
ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)
data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
ParEst_DM = t(apply(data_mat, 2, function(data){
  n_i = tabulate(data, nbins = M)
  idx_obs = which(n_i > 0)
  Kn = length(idx_obs)
  data_obs = n_i[idx_obs]
  
  #c) Dirichlet-Multinomial --> M = M
  # Param. estimation (DirMulti)
  start_params <- c(gamma = 0.5)
  fit <- optim(par = start_params, fn = llik_DirMult,
               n = n, M = M, data = n_i, # extra parameters
               method = "L-BFGS-B",
               lower = c(1e-10), upper = c(Inf))
  c(fit$par[1],fit$value)
})) 

round(ParEst_DM,0)
hh = 26
ParEst_DM[hh,]
n_i = tabulate(data_mat[,hh], M)
gamma_grid = seq(0.1, 1000, length.out = 5000)
llik_grid = rep(0,length(gamma_grid))
for(idx_gamma in seq_along(gamma_grid)){
  llik_grid[idx_gamma] = -llik_DirMult(gamma_grid[idx_gamma], n=n, M=M, data = n_i)
}

if(save_img)
  pdf("img/NegBin_llik_DM.pdf", width = width, height = height)
par( mfrow = c(1,1), mar = c(2.5,3,1,1), mgp=c(1.5,0.55,0), bty = "l", las = 1, cex = 2 )
plot(gamma_grid, llik_grid, type="l", lwd=3,
     ylim = quantile(llik_grid, c(0.02,1)),
     xlab="", ylab="")
mtext(expression(gamma), side=1, line=1.25, cex = 2)  # control xlab distance here
abline(h = -ParEst_DM[hh,2], lty = 2, lwd = 3, col = "red")
if(save_img)
  dev.off()