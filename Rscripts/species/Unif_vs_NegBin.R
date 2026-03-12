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
# Uniforme ----------------------------------------------------------------
name = "Uniform"
params = NA
ptrue = sim_generic_species(name,M,params)
ptrue = sort(ptrue, decreasing = TRUE)
ptrue_mat = matrix(nrow = Nrep,ncol = M)
ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)
data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
Nj_mat = apply(data_mat, 2, tabulate, nbins = M) # M x Nrep
Mmax_plot = M
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

EB_Unif = t(apply(data_mat, 2, function(data){
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
ExpMstar = apply(EB_Unif, 1, function(x){
  gamma = x[4]; Lambda = x[3]; Kn = x[6]
  Mstar_ub = 5000
  ## Mstar
  logExpMstar = compute_logV( Kn+1, n, gamma, Lambda, M_max ) - 
    compute_logV( Kn,   n, gamma, Lambda, M_max )
  exp(logExpMstar) # Expected value
})
EB_Unif = cbind(EB_Unif, ExpMstar)
colnames(EB_Unif) = c("sigma","theta","LambdaFDP","gammaFDP","gammaDM","Kn",
                      "Mmax","Freq",
                      "UB_PD","UB_FDP","UB_DM","pmin_obs","Mstar")

# NegBin ----------------------------------------------------------------
name = "NegBin"
params = c(1,0.003)
Nrep = 100
ptrue = sim_generic_species(name,M,params)
ptrue = sort(ptrue, decreasing = TRUE)
ptrue_mat = matrix(nrow = Nrep,ncol = M)
ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)
data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
Nj_mat = apply(data_mat, 2, tabulate, nbins = M) # M x Nrep
Mmax_plot = M
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


EB_NegBin = t(apply(data_mat, 2, function(data){
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
ExpMstar = apply(EB_NegBin, 1, function(x){
  gamma = x[4]; Lambda = x[3]; Kn = x[6]
  Mstar_ub = 5000
  ## Mstar
  logExpMstar = compute_logV( Kn+1, n, gamma, Lambda, M_max ) - 
    compute_logV( Kn,   n, gamma, Lambda, M_max )
  exp(logExpMstar) # Expected value
})
EB_NegBin = cbind(EB_NegBin, ExpMstar)
colnames(EB_NegBin) = c("sigma","theta","LambdaFDP","gammaFDP","gammaDM","Kn",
                      "Mmax","Freq",
                      "UB_PD","UB_FDP","UB_DM","pmin_obs","Mstar")

# Comparison --------------------------------------------------------------
EB_Unif[,7:12] = EB_Unif[,7:12]*1000
EB_NegBin[,7:12] = EB_NegBin[,7:12]*1000
round(colMeans(EB_Unif),1); round(colMeans(EB_NegBin),1)
round(apply(EB_Unif, 2, sd),2); round(apply(EB_NegBin, 2, sd),2)



EB_NegBin[,12] <= EB_NegBin[,7]

# Coverage ----------------------------------------------------------------

Cov_Unif = t(apply(EB_Unif, 1, function(x){
  Mmax = x[7]
  UBs = x[8:11]
  Mmax <= UBs
}))
colSums(Cov_Unif)/Nrep

Cov_NegBin = t(apply(EB_NegBin, 1, function(x){
  Mmax = x[7]
  UBs = x[8:11]
  Mmax <= UBs
}))
colSums(Cov_NegBin)/Nrep




# Plots -------------------------------------------------------------------


# Est. Quantities ---------------------------------------------------------
ymax_all = c(215,20000,20000,230)
ymin_all = c(190,0,0,180)
idx = 1
for(v in c(3,4,5,13)){
  xx = EB_Unif[,v]; yy = EB_NegBin[,v]
  if(v == 13){
    xx = xx + EB_Unif[,6]
    yy = yy + EB_NegBin[,6]
  }
  xpos = c(1,5); xlabs = c("Unif","NegBin")
  ymax = ymax_all[idx];ymin = ymin_all[idx]
  ylim_plot = c(ymin,ymax)
  ylab = colnames(EB_NegBin)[v]
  
  par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
  plot(1,1,type = "n",xlim = c(0.5,6), ylim = ylim_plot,
       xlab = "r", ylab = ylab, 
       xaxt = "n", yaxt = "n",
       main = "",
       cex.lab = cex.lab, cex.axis = cex.axis)
  axis(side = 2, las = 1, cex.axis = cex.axis)
  axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
  grid(lty = 1,lwd = 1, col = "gray90" )
  boxplot(xx, at = xpos[1], add = T, 
          col = "pink", pch = 16, yaxt = "n", cex = 0.5)
  boxplot(yy, at = xpos[2], add = T, 
          col = "skyblue", pch = 16, yaxt = "n", cex = 0.5)
  if(v == 13)
    abline(h = M, lty = 3, col = "red", lwd = 2)
  
  idx = idx + 1
}

# UB -------------------------------------------------------------------
ymax = c(12)
ymin = 4.5
ylim_plot = c(ymin,ymax)
ylab = "1000*bound"
pos_unif = c(0.5,1,1.5)
pos_nb   = c(4.5,5,5.5)
xlabs = c("Unif","NegBin"); xpos = c(1,5)
  
par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(0,6), ylim = ylim_plot,
     xlab = "r", ylab = ylab, 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
grid(lty = 1,lwd = 1, col = "gray90" )
idx = 1
for(v in c(7,10,11)){
  xx = EB_Unif[,v]; yy = EB_NegBin[,v]
  boxplot(xx, at = pos_unif[idx], add = T, 
          col = "pink", pch = 16, yaxt = "n", cex = 0.5)
  boxplot(yy, at = pos_nb[idx], add = T, 
          col = "skyblue", pch = 16, yaxt = "n", cex = 0.5)
  
  
  idx = idx + 1
}



# UB - function of params -------------------------------------------------------------------
UB_DM_fix = rep(1,Nrep)
gamma = 100
for(ii in 1:Nrep){
  Kn = EB_NegBin[ii,6]
  ubDM = exp(compute_log_UB_DirMulti( Rmax, gamma, M, Kn, n, alpha ))
  ubDM = min(ubDM,1)
  UB_DM_fix[ ii] = ubDM
}
# --> se gamma enorme, va a 1/M perché vuol dire che sta imponendo che 
# sia il caso uniforme
# --> Se gamma è troppo grande non copre perché in realtà le prob non sono
# uniformi, Mmax è spesso maggiore di 1/M


# Caso FDP
UB_FDP_fix = rep(1,Nrep)
Lambda = 50
for(ii in 1:Nrep){
  Kn = EB_NegBin[ii,6]
  ubFD = exp(compute_log_UBMarkov_FD( Rmax, gamma, Lambda, Kn, n, alpha, M_max ))
  ubFD = min(ubFD,1)
  UB_FDP_fix[ ii] = ubFD
}


Nrep = 1

ymax = c(12)
ymin = 4.5
ylim_plot = c(ymin,ymax)
ylab = "1000*bound"
pos_unif = c(0.5,1,1.5)
pos_nb   = c(4.5,5,5.5,4,6)
xlabs = c("Unif","NegBin"); xpos = c(1,5)

par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(0,7), ylim = ylim_plot,
     xlab = "r", ylab = ylab, 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
grid(lty = 1,lwd = 1, col = "gray90" )
idx = 1
for(v in c(7,10,11)){
  xx = EB_Unif[,v]; yy = EB_NegBin[,v]
  boxplot(xx, at = pos_unif[idx], add = T, 
          col = "pink", pch = 16, yaxt = "n", cex = 0.5)
  boxplot(yy, at = pos_nb[idx], add = T, 
          col = "skyblue", pch = 16, yaxt = "n", cex = 0.5)
  
  
  idx = idx + 1
}

boxplot(UB_FDP_fix*1000, at = pos_nb[idx], add = T, 
        col = "red", pch = 16, yaxt = "n", cex = 0.5)
idx = idx + 1
boxplot(UB_DM_fix*1000, at = pos_nb[idx], add = T, 
        col = "blue", pch = 16, yaxt = "n", cex = 0.5)



# gamma varies - Lambda grid
Lambda_grid = c(50,150,200,250,300,500)
gamma_grid = c(seq(0.001,1,length.out = 100),
               seq(1,150,length.out = 200))
UB_FDP_list = vector("list", length = length(gamma_grid))
res = vector("list", length = length(Lambda_grid))
Nrep = 1
counter = 1
for(gg in seq_along(Lambda_grid)){
  Lambda = Lambda_grid[gg]
  for(hh in seq_along(gamma_grid)) {
    UB_FDP_fix = rep(1,Nrep)
    gamma = gamma_grid[hh]
    for(ii in 1:Nrep){
      Kn = EB_NegBin[ii,6]
      ubFDP = exp(compute_log_UBMarkov_FD( Rmax, gamma, Lambda, Kn, n, alpha, M_max ))
      ubFDP = min(ubFDP,1)
      UB_FDP_fix[ii] = ubFDP
    }
    UB_FDP_list[[hh]] = UB_FDP_fix
    counter = counter + 1
  }
  res[[gg]] = 1000*sapply(UB_FDP_list, mean)
}

par( mfrow = c(2,3), mar = c(3.5,2.5,1,1), mgp=c(1,1,0), bty = "l", las = 1, cex.lab = cex.lab )
for(gg in seq_along(Lambda_grid)){
  plot(x = gamma_grid, y = res[[gg]], type = "b", 
       pch = 16, ylab = Lambda_grid[gg], ylim = c(2,12.5))
  abline(h = 1000/M, lty = 4, col = "red")
  abline(h = EB_NegBin[,7], lty = 4, col = "green")
  # abline(v = EB_NegBin[,5], lty = 1, lwd = 1, col = "skyblue")
}


# Worst-Uniform case ----------------------------------------------------
n = 500 
name = "WorstUnif"
params = c(n,alpha)
ptrue = sim_generic_species(name,M,params)
ptrue = sort(ptrue, decreasing = TRUE)
ptrue_mat = matrix(nrow = Nrep,ncol = M)
ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)
data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
Nj_mat = apply(data_mat, 2, tabulate, nbins = M) # M x Nrep
Mmax_plot = M
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

EB_WorstUnif = t(apply(data_mat, 2, function(data){
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
ExpMstar = apply(EB_WorstUnif, 1, function(x){
  gamma = x[4]; Lambda = x[3]; Kn = x[6]
  Mstar_ub = 5000
  ## Mstar
  logExpMstar = compute_logV( Kn+1, n, gamma, Lambda, M_max ) - 
    compute_logV( Kn,   n, gamma, Lambda, M_max )
  exp(logExpMstar) # Expected value
})
EB_WorstUnif = cbind(EB_WorstUnif, ExpMstar)
colnames(EB_WorstUnif) = c("sigma","theta","LambdaFDP","gammaFDP","gammaDM","Kn",
                      "Mmax","Freq",
                      "UB_PD","UB_FDP","UB_DM","pmin_obs","Mstar")



EB_WorstUnif[,7:12] = EB_WorstUnif[,7:12]*1000
round(colMeans(EB_WorstUnif),1)
Cov_WorstUnif = t(apply(EB_WorstUnif, 1, function(x){
  Mmax = x[7]
  UBs = x[8:11]
  Mmax <= UBs
}))
colSums(Cov_WorstUnif)/Nrep



ymax_all = c(75,200,1,75)
ymin_all = c(65,0,0,65)
idx = 1
for(v in c(3,4,5,13)){
  xx = EB_WorstUnif[,v]
  if(v == 13){
    xx = xx + EB_WorstUnif[,6]
  }
  xpos = c(3); xlabs = c("WorstUnif")
  ymax = ymax_all[idx];ymin = ymin_all[idx]
  ylim_plot = c(ymin,ymax)
  ylab = colnames(EB_WorstUnif)[v]
  
  par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
  plot(0,0,type = "n",xlim = c(0,6), ylim = ylim_plot,
       xlab = "r", ylab = ylab, 
       xaxt = "n", yaxt = "n",
       main = "",
       cex.lab = cex.lab, cex.axis = cex.axis)
  axis(side = 2, las = 1, cex.axis = cex.axis)
  axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
  grid(lty = 1,lwd = 1, col = "gray90" )
  boxplot(xx, at = xpos[1], add = T, 
          col = "green", pch = 16, yaxt = "n", cex = 0.5)
  if(v == 13)
    abline(h = M, lty = 3, col = "red", lwd = 2)
  
  idx = idx + 1
}



round(EB_WorstUnif[,c(3,4,5,6,7,10,11,13)],2)
summary(EB_WorstUnif[,c(3,4)])


# --> DirMulti fa male perché dico M=200 quando M=70. Quindi gamma è molto piccolo
# per riuscire a spiegare come mai le prob. vanno a zero ma questo fa saltare tutto
# --> FDP ci prende perché riesce e stimare correttamente la dimensione dell'alfabeto
# e poi sappiamo che con le prob uniformi va abbastanza bene

UB_DM_fix = rep(1,Nrep)
gamma = 10
for(ii in 1:Nrep){
  Kn = EB_WorstUnif[ii,6]
  ubDM = exp(compute_log_UB_DirMulti( Rmax, gamma, M, Kn, n, alpha ))
  ubDM = min(ubDM,1)
  UB_DM_fix[ ii] = ubDM
}
UB_DM_fix = UB_DM_fix * 1000

UB_FDP_fix = rep(1,Nrep)
Lambda = 70
for(ii in 1:Nrep){
  Kn = EB_WorstUnif[ii,6]
  ubFD = exp(compute_log_UBMarkov_FD( Rmax, gamma, Lambda, Kn, n, alpha, M_max ))
  ubFD = min(ubFD,1)
  UB_FDP_fix[ ii] = ubFD
}
UB_FDP_fix = UB_FDP_fix*1000

1000/70
summary(UB_DM_fix)
summary(UB_FDP_fix)
summary(EB_WorstUnif[,10])
summary(EB_WorstUnif[,11])




# NegBin - DM fails ------------------------------------------------------------------


gamma_grid = c(seq(0.001,1,length.out = 500),
               seq(1,3000,length.out = 1000))
UB_DM_list = vector("list", length = length(gamma_grid))
for(hh in seq_along(gamma_grid)) {
  UB_DM_fix = rep(1,Nrep)
  gamma = gamma_grid[hh]
  for(ii in 1:Nrep){
    Kn = EB_NegBin[ii,6]
    ubDM = exp(compute_log_UB_DirMulti( Rmax, gamma, M, Kn, n, alpha ))
    ubDM = min(ubDM,1)
    UB_DM_fix[ii] = ubDM
  }
  UB_DM_list[[hh]] = UB_DM_fix
}

res = 1000*sapply(UB_DM_list, mean)
plot(x = gamma_grid, y = res, type = "b", pch = 16)
abline(h = 1000/M, lty = 4, col = "red")
abline(h = EB_NegBin[,7], lty = 4, col = "green")
abline(v = EB_NegBin[,5], lty = 1, lwd = 1, col = "skyblue")


# Unif - DM Ok ------------------------------------------------------------------


gamma_grid = c(seq(0.001,1,length.out = 500),
               seq(1,3000,length.out = 1000))
UB_DM_list = vector("list", length = length(gamma_grid))
for(hh in seq_along(gamma_grid)) {
  UB_DM_fix = rep(1,Nrep)
  gamma = gamma_grid[hh]
  for(ii in 1:Nrep){
    Kn = EB_Unif[ii,6]
    ubDM = exp(compute_log_UB_DirMulti( Rmax, gamma, M, Kn, n, alpha ))
    ubDM = min(ubDM,1)
    UB_DM_fix[ii] = ubDM
  }
  UB_DM_list[[hh]] = UB_DM_fix
}

res = 1000*sapply(UB_DM_list, mean)
plot(x = gamma_grid, y = res, type = "b", pch = 16)
abline(h = 1000/M, lty = 4, col = "red")
abline(h = EB_Unif[,7], lty = 4, col = "green")
abline(v = EB_Unif[,5], lty = 1, lwd = 1, col = "skyblue")



# WorstUnif - DM fails ------------------------------------------------------------------


gamma_grid = c(seq(0.001,1,length.out = 500),
               seq(1,10,length.out = 500))
UB_DM_list = vector("list", length = length(gamma_grid))
for(hh in seq_along(gamma_grid)) {
  UB_DM_fix = rep(1,Nrep)
  gamma = gamma_grid[hh]
  for(ii in 1:Nrep){
    Kn = EB_WorstUnif[ii,6]
    ubDM = exp(compute_log_UB_DirMulti( Rmax, gamma, M, Kn, n, alpha ))
    ubDM = min(ubDM,1)
    UB_DM_fix[ii] = ubDM
  }
  UB_DM_list[[hh]] = UB_DM_fix
}

res = 1000*sapply(UB_DM_list, mean)
plot(x = gamma_grid, y = res, type = "b", pch = 16, ylim = c(3.5,15))
abline(h = 1000/M, lty = 4, col = "red")
abline(h = EB_WorstUnif[,7], lty = 4, col = "green")
abline(v = EB_WorstUnif[,5], lty = 1, lwd = 1, col = "skyblue")

