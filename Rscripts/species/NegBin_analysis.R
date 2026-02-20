# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
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
mycol = c("darkorange","darkred","darkblue","lightpink","aquamarine")
mycol2 = c("black","lightblue")
lgd_names = c("Freq","PD","FDP","DirMulti-m","DirMulti")
tcol = 0.25
# Options ---------------------------------------------------------------------
params_zipfs = list(0.9,1.02,2)
params_geom = list(0.85,0.9,0.95)
params_unif = list(NA)
params_negbin = list(c(1,0.003),c(5,0.003),c(1,0.01))
experiments = list("Zipfs" = params_zipfs,
                   "Geom" = params_geom,
                   "Uniform" = params_unif,
                   "NegBin" = params_negbin)
Nrep = 500
n = 500
M = 200 
M_max = 200
seed = 232131
set.seed(seed)
Mmax = M;
AccCrv_length = 50
ii = 4; name = names(experiments)[[ii]]
jj = 1
params = experiments[[ii]][[jj]]
trim_params = sapply(params, get_first3digits, 4)
if(length(trim_params)>1)
  trim_params = paste0(trim_params[1],"_",trim_params[2])


save_name_base = paste0("save/Species_") 
save_name_base_cov = paste0("save/Species_Cov_") 
img_fld = paste0("img/Species_") 

filename = paste0(save_name_base,name,"_nfix_",trim_params,".Rdat")
load(filename)
# Missing mass with M -----------------------------------------------------

Mmin_grid = 50; Mmax_grid = 1000
Mgrid = seq(Mmin_grid,Mmax_grid,by = 50); LMgrid = length(Mgrid)
Nexp = length(Mgrid)

MissMass_list = vector("list",LMgrid)
for(ii in 1:LMgrid){
  M = Mgrid[ii]
  ptrue = sim_generic_species(name,M,params)
  ptrue = sort(ptrue, decreasing = TRUE)
  ptrue_mat = matrix(nrow = Nrep,ncol = M)
  ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)
  
  data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
  Missing_mass = apply(rbind(data_mat,ptrue_mat), 2, function(x){
    data = x[1:n]
    pj = x[(n+1):(n+M)]
    n_i = tabulate(data, nbins = M)
    idx_no_obs = which( n_i == 0 )
    sum(pj[idx_no_obs])
  }) # M x Nrep
  
  MissMass_list[[ii]] = Missing_mass
}



pos = 1:LMgrid
xpos = pos; xlabs = Mgrid
ymax = max(sapply(MissMass_list,max));ymin = min(sapply(MissMass_list,min))
ylim_plot = c(ymin,ymax)
ylab = "Missing mass"

img_name = paste0("img/",name,"_",trim_params,"_MissingMass",".pdf")
if(save_img)
  pdf(img_name, width = width+6, height = height)
par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(0,(LMgrid+1)), ylim = ylim_plot,
       xlab = "r", ylab = ylab, 
       xaxt = "n", yaxt = "n",
       main = "",
       cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
grid(lty = 1,lwd = 1, col = "gray90" )
for(ii in 1:(LMgrid)){
  boxplot(MissMass_list[[ii]], at = pos[ii], add = T, 
          col = "grey90", pch = 16, yaxt = "n", cex = 0.5)
}
if(save_img)
  dev.off()



# Mmax with M -------------------------------------------------------------
pos = 1:LMgrid
xpos = pos; xlabs = Mgrid
ymax = max(sapply(ExpRes_list,function(x) max(x[,1])));ymin = min(sapply(ExpRes_list,function(x) min(x[,1])))
ylim_plot = c(ymin,ymax)
ypos = seq(ymin,ymax,length.out = 5); ylabs = round(ypos*1e3,0)
ylab = "M_max"

img_name = paste0("img/",name,"_",trim_params,"_Mmax",".pdf")
if(save_img)
  pdf(img_name, width = width+6, height = height)
par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(0,(LMgrid+1)), ylim = ylim_plot,
     xlab = "r", ylab = ylab, 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
grid(lty = 1,lwd = 1, col = "gray90" )
for(ii in 1:(LMgrid)){
  boxplot(ExpRes_list[[ii]][,1], at = pos[ii], add = T, 
          col = "grey90", pch = 16, yaxt = "n", cex = 0.5)
}
if(save_img)
  dev.off()


# Estimates - EB ---------------------------------------------------------------
Mmin_grid = 50; Mmax_grid = 1000
Mgrid = seq(Mmin_grid,Mmax_grid,by = 50); LMgrid = length(Mgrid)
Nexp = length(Mgrid)
Nrep = 100
ParEst_list = vector("list",LMgrid)
for(ii in 1:LMgrid){
  cat("\n",ii,"/",LMgrid,"\n")
  M = Mgrid[ii]
  ptrue = sim_generic_species(name,M,params)
  ptrue = sort(ptrue, decreasing = TRUE)
  ptrue_mat = matrix(nrow = Nrep,ncol = M)
  ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)
  
  data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
  ParEst = t(apply(data_mat, 2, function(data){
    n_i = tabulate(data, nbins = M)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    data_obs = n_i[idx_obs]
    
    res = matrix(nrow = 1, ncol = 5)
    colnames(res) = c("sigma-PD","theta-PD","Lambda-FDP","gamma-FDP","gamma-DirMulti")
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
    
    res
  })) # Nrep x 3 matrix
  ParEst_list[[ii]] = ParEst
}

gammas_DM = lapply(ParEst_list, function(ParEst) ParEst[,5])

pos = 1:LMgrid
xpos = pos; xlabs = Mgrid
ymax = 1000;ymin = 0
ylim_plot = c(ymin,ymax)
ypos = seq(ymin,ymax,length.out = 5); ylabs = round(ypos,2)
ylab = "Gamma (DirMulti)"

img_name = paste0("img/",name,"_",trim_params,"_Mmax",".pdf")
if(save_img)
  pdf(img_name, width = width+6, height = height)
par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(0,(LMgrid+1)), ylim = ylim_plot,
     xlab = "r", ylab = ylab, 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
grid(lty = 1,lwd = 1, col = "gray90" )
for(ii in 1:(LMgrid)){
  boxplot(gammas_DM[[ii]], at = pos[ii], add = T, 
          col = "grey90", pch = 16, yaxt = "n", cex = 0.5)
}
if(save_img)
  dev.off()


# Estimates - MCMC ---------------------------------------------------------------
Mmin_grid = 50; Mmax_grid = 500
Mgrid = seq(Mmin_grid,Mmax_grid,by = 50); LMgrid = length(Mgrid)
Nexp = length(Mgrid)
Nrep = 10
ParEst_list = vector("list",LMgrid)
for(ii in 1:LMgrid){
  cat("\n",ii,"/",LMgrid,"\n")
  M = Mgrid[ii]
  ptrue = sim_generic_species(name,M,params)
  ptrue = sort(ptrue, decreasing = TRUE)
  ptrue_mat = matrix(nrow = Nrep,ncol = M)
  ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)
  
  data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
  ParEst = t(apply(data_mat, 2, function(data){
    n_i = tabulate(data, nbins = M)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    data_obs = n_i[idx_obs]
    
    res = matrix(nrow = 1, ncol = 5)
    colnames(res) = c("sigma-PD","theta-PD","Lambda-FDP","gamma-FDP","gamma-DirMulti")
    
    model = "DirMulti"
    Niter = 10000
    init_val = c(1)
    hy_prior = c(1,1) 
    Adp_var = c(0.1)
    UpdateParam = c(TRUE)
    fit = ParEst_MCMC_generic(model=model,n=n,Kn=Kn,Nj=n_i,Niter=Niter,
                              init_val=init_val,hy_prior=hy_prior,Adp_var=Adp_var,UpdateParam=UpdateParam,
                              M=M,seed = seed, M_max = 500)
    gammas = fit$gamma_mcmc
    res[1,5] = mean(gammas[(Niter/2):Niter])    
    res
  })) # Nrep x 3 matrix
  ParEst_list[[ii]] = ParEst
}

gammas_DM = lapply(ParEst_list, function(ParEst) ParEst[,5])

pos = 1:LMgrid
xpos = pos; xlabs = Mgrid
ymax = 10;ymin = 0
ylim_plot = c(ymin,ymax)
ypos = seq(ymin,ymax,length.out = 5); ylabs = round(ypos,2)
ylab = "Gamma (DirMulti)"

img_name = paste0("img/",name,"_",trim_params,"_Mmax",".pdf")
if(save_img)
  pdf(img_name, width = width+6, height = height)
par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(0,(LMgrid+1)), ylim = ylim_plot,
     xlab = "r", ylab = ylab, 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
grid(lty = 1,lwd = 1, col = "gray90" )
for(ii in 1:(LMgrid)){
  boxplot(gammas_DM[[ii]], at = pos[ii], add = T, 
          col = "grey90", pch = 16, yaxt = "n", cex = 0.5)
}
if(save_img)
  dev.off()


# UB MCMC -----------------------------------------------------------------


Rmax = 100
alpha = 0.05
Mmin_grid = 50; Mmax_grid = 300
Mgrid = seq(Mmin_grid,Mmax_grid,by = 50); LMgrid = length(Mgrid)
Nexp = length(Mgrid)
Nrep = 100
UB_list = vector("list",LMgrid)
for(ii in 1:LMgrid){
  cat("\n",ii,"/",LMgrid,"\n")
  M = Mgrid[ii]
  ptrue = sim_generic_species(name,M,params)
  ptrue = sort(ptrue, decreasing = TRUE)
  ptrue_mat = matrix(nrow = Nrep,ncol = M)
  ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)
  
  data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
  UBres = t(apply(data_mat, 2, function(data){
    n_i = tabulate(data, nbins = M)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    data_obs = n_i[idx_obs]
    
    res = matrix(nrow = 1, ncol = 5)
    colnames(res) = c("sigma-PD","theta-PD","Lambda-FDP","gamma-FDP","gamma-DirMulti")
    
    model = "DirMulti"
    Niter = 10000
    init_val = c(1)
    hy_prior = c(1,1) 
    Adp_var = c(0.1)
    UpdateParam = c(TRUE)
    fit = ParEst_MCMC_generic(model=model,n=n,Kn=Kn,Nj=n_i,Niter=Niter,
                              init_val=init_val,hy_prior=hy_prior,Adp_var=Adp_var,UpdateParam=UpdateParam,
                              M=M,seed = seed, M_max = 500)
    gammas = fit$gamma_mcmc
    gamma = mean(gammas[(Niter/2):Niter])
    y = exp(compute_log_UB_DirMulti( Rmax, gamma, M, Kn, n, alpha ))
    y = min(y,1); y
    res[1,5] = y
    res
  })) # Nrep x 3 matrix
  UB_list[[ii]] = UBres
}

ub_DM_all = lapply(UB_list, function(UBres) UBres[,5])
oracle = sapply(ExpRes_list, function(x) quantile(x[,1], 1-alpha))


pos = 1:LMgrid
xpos = pos; xlabs = Mgrid
ymax = 14*1e-3;ymin = 0
ylim_plot = c(ymin,ymax)
ypos = seq(ymin,ymax,length.out = 5); ylabs = round(ypos*1e3,0)
ylab = "UB (DirMulti)"

img_name = paste0("img/",name,"_",trim_params,"_UBMCMC",".pdf")
if(save_img)
  pdf(img_name, width = width+6, height = height)
par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(0,(LMgrid+1)), ylim = ylim_plot,
     xlab = "r", ylab = ylab, 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
grid(lty = 1,lwd = 1, col = "gray90" )
points( x = pos, y = oracle[1:LMgrid], 
        type = "l", lwd = 5, col = "black" )
for(ii in 1:(LMgrid)){
  boxplot(ub_DM_all[[ii]], at = pos[ii], add = T, 
          col = "grey90", pch = 16, yaxt = "n", cex = 0.5)
}
if(save_img)
  dev.off()


sapply(1:LMgrid, function(ii){
  length(which(ub_DM_all[[ii]] < oracle[ii]))/length(ub_DM_all[[ii]])
})

ii = 2
dens_ub = density(ub_DM_all[[ii]])
plot(dens_ub$x,dens_ub$y, xlab = "", ylab = "")
