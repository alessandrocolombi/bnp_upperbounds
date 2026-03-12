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
width = 12; height = 8
cex.labels <- cex.lab <- 1.5
cex.axis <- 1.5
cex.legend <- 1.5
mycol = c("darkorange","darkred","darkblue","lightblue")
mycol_var = hcl.colors(n = 4, palette = "Greens", rev = TRUE)
mycol2 = c("black","lightblue")
lgd_names = c("Oracle","Freq","PD","FDP","Dir-Multi")


# Options -----------------------------------------------------------------
M = 200
n = 500
Rmax = 100
Nrep = 200
alpha <- alfa <- 0.05
seed = 42
M_max = 200

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

par( mfrow = c(1,1), mar = c(3.5,3.5,1,1), mgp=c(2,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(x = gamma_grid, y = U_DM_grid, type = "b", 
     pch = 16, ylab = "DirMulti", ylim = c(2,15), xlab = "gamma")


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

par( mfrow = c(2,3), mar = c(3.5,3.5,1,1), mgp=c(2,1,0), bty = "l", las = 1, cex.lab = cex.lab )
for(gg in seq_along(Lambda_grid)){
  plot(x = gamma_grid, y = res[[gg]], type = "b", 
       pch = 16, ylab = paste0("FDP - ",Lambda_grid[gg]), ylim = c(2,12.5), 
       xlab = "gamma")
}


## DP - bivariate -------------------------------------------------------
Kn = 70 # -> la forma non cambia ma il massimo si, in funzione di Kn
sigma_grid = c(0,0.25,0.5,0.7,0.9,0.99)
theta_grid = c(seq(0.001,1,length.out = 100),
               seq(1,150,length.out = 200))
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

par( mfrow = c(2,3), mar = c(3.5,3.5,1,1), mgp=c(2,1,0), bty = "l", las = 1, cex.lab = cex.lab )
for(gg in seq_along(sigma_grid)){
  plot(x = theta_grid, y = res[[gg]], type = "b", 
       pch = 16, ylab = paste0("PD - ",sigma_grid[gg]), ylim = c(2,12.5), 
       xlab = "theta")
}
