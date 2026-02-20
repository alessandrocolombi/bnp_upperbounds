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

# # From BinomialCIs
# source("../../../BinomialCIs/R/Rfunctions.R")
# Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")



# Custom functions --------------------------------------------------------

# Test DirMulti -----------------------------------------------------------
model = "DirMulti"
M <- Mmax <- 200
n = 100
Nrep = 5
gamma_DM = 0.5
params = c(gamma_DM,M)
seed = 312312
set.seed(seed)
ptrue_mat = SimModel_generic(model,Mmax,Nrep,params,seed)

plot(ptrue_mat[,1], pch = 16, cex = 0.5, type = "p")

data_mat = SimData(n,ptrue_mat,seed)

data = data_mat[,1]
n_i = tabulate(data, nbins = M)
idx_obs = which(n_i > 0)
Kn = length(idx_obs)
Nj_obs = n_i[idx_obs]

model = "DirMulti"
Niter = 5000
init_val = c(1)
hy_prior = c(1,0.01) 
Adp_var = c(0.1)
UpdateParam = c(TRUE)
fit = ParEst_MCMC_generic(model=model,n=n,Kn=Kn,Nj=n_i,Niter=Niter,
                          init_val=init_val,hy_prior=hy_prior,Adp_var=Adp_var,UpdateParam=UpdateParam,
                          M=M,seed = seed, M_max = 500)
gammas = fit$gamma_mcmc
plot(gammas, type = "l", lwd = 0.5)
mean(gammas)
median(gammas)


# Test FDP -----------------------------------------------------------
model = "FDP"
M <- Mmax <- 200
n = 1000
Nrep = 5
gamma_FDP = 0.5
Lambda_FDP = 150
params = c(gamma_FDP,Lambda_FDP)
seed = 312312
set.seed(seed)
ptrue_mat = SimModel_generic(model,Mmax,Nrep,params,seed)

plot(ptrue_mat[,1], pch = 16, cex = 0.5, type = "p")

data_mat = SimData(n,ptrue_mat,seed)

data = data_mat[,1]
n_i = tabulate(data, nbins = M)
idx_obs = which(n_i > 0)
Kn = length(idx_obs)
Nj_obs = n_i[idx_obs]

model = "FDP"
Niter = 5000
init_val = c(1,100)
hy_prior = c(1,1,1,1) 
Adp_var = c(0.1,0.1)
UpdateParam = c(TRUE,TRUE)
fit = ParEst_MCMC_generic(model=model,n=n,Kn=Kn,Nj=Nj_obs,Niter=Niter,
                          init_val=init_val,hy_prior=hy_prior,Adp_var=Adp_var,UpdateParam=UpdateParam,
                          M=M,seed = seed, M_max = 500)
gammas = fit$gamma_mcmc
lambdas = fit$Lambda_mcmc
plot(gammas, type = "l", lwd = 0.5)
plot(lambdas, type = "l", lwd = 0.5)
mean(gammas)
median(gammas)

# PD ----------------------------------------------------------------------
model = "PD"
M <- Mmax <- 200
n = 1000
Nrep = 5
sigma_PD = 0.15
theta_PD = 10
params = c(sigma_PD,theta_PD)
seed = 312312
set.seed(seed)
ptrue_mat = SimModel_generic(model,Mmax,Nrep,params,seed)

plot(ptrue_mat[,1], pch = 16, cex = 0.5, type = "p")

data_mat = SimData(n,ptrue_mat,seed)

data = data_mat[,1]
n_i = tabulate(data, nbins = M)
idx_obs = which(n_i > 0)
Kn = length(idx_obs)
Nj_obs = n_i[idx_obs]

log_eppfPYP(n,Kn,Nj_obs,sigma_PD,theta_PD)
log_eppfPYP(n,Kn,Nj_obs,0.5,theta_PD)
log_eppfPYP(n,Kn,Nj_obs,0.75,theta_PD)

log_eppfPYP(n,Kn,Nj_obs,sigma_PD,theta_PD) > log_eppfPYP(n,Kn,Nj_obs,0.5,theta_PD)
log_eppfPYP(n,Kn,Nj_obs,sigma_PD,theta_PD) > log_eppfPYP(n,Kn,Nj_obs,0.75,theta_PD)

model = "PD"
Niter = 5000
init_val = c(0.5,20)
hy_prior = c(1,1,20,1) 
Adp_var = c(1,1)
UpdateParam = c(TRUE,TRUE)
fit = ParEst_MCMC_generic(model=model,n=n,Kn=Kn,Nj=Nj_obs,Niter=Niter,
                          init_val=init_val,hy_prior=hy_prior,Adp_var=Adp_var,UpdateParam=UpdateParam,
                          M=M,seed = seed, M_max = 500)

fit$sigma_mcmc
mean(fit$sigma_mcmc)
fit$theta_mcmc
mean(fit$theta_mcmc)


# Param. estimation (PYP)
start_params <- c(alpha = 0.5, theta = 1)
fit <- optim(par = start_params, fn = llik_pyp, 
             n = n, Kn = Kn, data_obs = Nj_obs, # extra parameters
             method = "L-BFGS-B",
             lower = c(0, -1), upper = c(1-1e-10, Inf)) 
alpha_mle = fit$par[1]
alpha_mle
theta_mle = fit$par[2]
theta_mle
fit$value
