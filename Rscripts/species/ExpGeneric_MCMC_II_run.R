# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[4] # <--- modify here
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
seed = 121321
M = 100
ExpGeneric_speciesMCMC_nfix_run = function(M, n, name, params, var_prior,
                                           alpha = 0.05,Rmax = 100, M_max = 200, 
                                           M_DM = NULL, seed = 121321)
{
  source("../../R/Rfunctions.R")
  Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
  # From BinomialCIs
  source("../../../BinomialCIs/R/Rfunctions.R")
  Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")
  
  set.seed(seed)
  ptrue = sim_generic_species(name,M,params)
  ptrue = sort(ptrue, decreasing = TRUE)
  # Define return object
  res_names = c("Mmax","Freq","PD","FDP","DirMulti",
                "sigma","theta",
                "gamma_FDP","Lambda_FDP",
                "gamma_DM", "Kn")
  res = matrix(NA,nrow = 1, ncol = length(res_names))
  colnames(res) = res_names
  #a) Generate data
  data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
  n_i = tabulate(data, nbins = M)
  idx_obs = which(n_i > 0)
  Kn = length(idx_obs)
  res[1,11] = Kn
  data_obs = n_i[idx_obs]
  if(Kn == M){
    Mmax = 0
  }else{
    idx_unobs = which(n_i == 0)
    Mmax = max(ptrue[idx_unobs])
  }
  #b) Freq
  pain = ub_pain(n = n, Rmax = Rmax, alfa = alpha)
  pain = min(pain,1)
  
  ### Prior definition - MCMC params
  Niter = 10000
  
  ## a) PD hyperparameters
  a_sigma = 1
  b_sigma = 1
  mu_theta = Hychoice_MCMC_general("PD", n, Kn, hy_prior = c(a_sigma,b_sigma) )
  a_theta = mu_theta*mu_theta/var_prior
  b_theta = mu_theta/var_prior
  
  ## b) FDP hyperparameters
  mu_Lambda = Kn 
  a_Lambda = mu_Lambda*mu_Lambda/var_prior
  b_Lambda = mu_Lambda/var_prior
  mu_gamma = 1 #Hychoice_MCMC_general("FDP", n, Kn, hy_prior = c(mu_Lambda) )
  a_gamma  = mu_gamma*mu_gamma/var_prior
  b_gamma  = mu_gamma/var_prior
  
  ## c) DirMulti hyperparameters
  mu_gammaDM = 1
  a_gammaDM  = mu_gammaDM*mu_gammaDM/var_prior
  b_gammaDM  = mu_gammaDM/var_prior
  
  ### Bayesian UB computation
  #c) Upper bound (PD)
  model = "PD"
  init_val = c(0.5,Kn)
  hy_prior = c(a_sigma,b_sigma,a_theta,b_theta) 
  Adp_var = c(0.1,0.1)
  UpdateParam = c(TRUE,TRUE)
  fit = ParEst_MCMC_generic(model=model,n=n,Kn=Kn,Nj=data_obs,Niter=Niter,
                            init_val=init_val,hy_prior=hy_prior,Adp_var=Adp_var,UpdateParam=UpdateParam,
                            M=M,seed = seed, M_max = 500)
  sigmas = fit$sigma_mcmc
  thetas = fit$theta_mcmc
  sigma  = mean(sigmas[(Niter/2):Niter])
  theta = mean(thetas[(Niter/2):Niter])
  ubpyp = exp(compute_log_UBMarkov( Rmax, sigma, theta, Kn, n, alpha ))
  ubpyp = min(ubpyp,1)
  res[1,6] = sigma
  res[1,7] = theta
  #c) Upper bound (FDP)
  model = "FDP"
  init_val = c(1,Kn)
  hy_prior = c(a_gamma,b_gamma,a_Lambda,b_Lambda) 
  Adp_var = c(0.1,0.1)
  UpdateParam = c(TRUE,TRUE)
  fit = ParEst_MCMC_generic(model=model,n=n,Kn=Kn,Nj=data_obs,Niter=Niter,
                            init_val=init_val,hy_prior=hy_prior,Adp_var=Adp_var,UpdateParam=UpdateParam,
                            M=M,seed = seed, M_max = 500)
  gammas = fit$gamma_mcmc
  Lambdas = fit$Lambda_mcmc
  gamma = mean(gammas[(Niter/2):Niter])
  Lambda = mean(Lambdas[(Niter/2):Niter])
  res[1,8] = gamma
  res[1,9] = Lambda
  ubFDP = exp(compute_log_UBMarkov_FD( Rmax, gamma, Lambda, Kn, n, alpha, M_max ))
  ubFDP = min(ubFDP,1)
  #d) Upper bound (DirMulti)
  model = "DirMulti"
  init_val = c(1)
  hy_prior = c(a_gammaDM,b_gammaDM) 
  Adp_var = c(0.1)
  UpdateParam = c(TRUE)
  fit = ParEst_MCMC_generic(model=model,n=n,Kn=Kn,Nj=n_i,Niter=Niter,
                            init_val=init_val,hy_prior=hy_prior,Adp_var=Adp_var,UpdateParam=UpdateParam,
                            M=M,seed = seed, M_max = 500)
  gammas = fit$gamma_mcmc
  gamma = mean(gammas[(Niter/2):Niter])
  res[1,10] = gamma
  ubDM = exp(compute_log_UB_DirMulti( Rmax, gamma, M, Kn, n, alpha ))
  ubDM = min(ubDM,1)
  #f) Save results
  res[1,1] = Mmax
  res[1,2] = pain
  res[1,3] = ubpyp
  res[1,4] = ubFDP
  res[1,5] = ubDM
  # Return
  res
}

ExpGeneric_speciesMCMC_nfix = function(name,M,n,params,Nrep = 100, 
                                       alpha = 0.05,Rmax = 100, 
                                       M_max = 200, seed0 = 121321,
                                       M_DM = NULL, var_prior = 1,
                                       parallel = TRUE, num_cores = 5)
{
  cat("\n M = ",M," || var_prior = ",var_prior,"\n")
  set.seed(seed0)
  seeds = sample(1:999999, size = Nrep)
  
  # Sequential case
  if(!parallel){
    res_list = lapply(seeds, function(seed) ExpGeneric_speciesMCMC_nfix_run(
      M=M, n=n, params=params, 
      alpha=alpha, var_prior=var_prior,
      Rmax=Rmax, M_max=M_max,
      M_DM=M_DM,
      seed=seed, name=name) )
  }else{
    ## Parallel case
    cluster <- makeCluster(num_cores, type = "SOCK")
    doSNOW::registerDoSNOW(cluster)
    clusterExport(cluster, list("alpha"), envir = environment())
    res_list = parLapply( cl = cluster, x = seeds,
                          fun = ExpGeneric_speciesMCMC_nfix_run,
                          M=M, n=n, alpha=alpha, name=name, params=params, 
                          M_DM=M_DM, var_prior=var_prior,
                          Rmax=Rmax, M_max=M_max )
    stopCluster(cluster)
  }
  
  res_mat = do.call(rbind,res_list)
  return(res_mat)
}

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
params_zipfs = list(0.9,1.02,2)
params_geom = list(0.85,0.9,0.95)
params_unif = list(NA)
params_negbin = list(c(1,0.003)) #list(c(1,0.003),c(5,0.003),c(1,0.01))
experiments = list("Zipfs" = params_zipfs,
                   "Geom" = params_geom,
                   "Uniform" = params_unif,
                   "NegBin" = params_negbin)

alpha <- alfa <- 0.05
num_cores = 33 # <---
Nrep = 200 # <---
n = 500
Rmax = 100; RmaxFD = 50
Mmin_grid = 50; Mmax_grid = 1000
Mgrid = seq(Mmin_grid,Mmax_grid,by = 50); LMgrid = length(Mgrid)
Nexp = length(Mgrid)
M_max = 200

seed0 = 42
set.seed(seed0)
seeds = sample(1:999999, size = Nexp)

save_exp = TRUE # <---
save_name_base = paste0("save/Species_MCMC_II_") 
img_fld = paste0("img/") 
# n fix -----------------------------------------------------------------
igrid = c(3,4)
ii = 4
for(ii in igrid){
  
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 1
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = sapply(params, get_first3digits, 4)
    if(length(trim_params)>1)
      trim_params = paste0(trim_params[1],"_",trim_params[2])
    
    run_obj <- Map(function(m, s) list(M = m, seed = s),Mgrid, seeds)
    
    ExpRes_list = lapply( run_obj,
                          function(x){
                            ExpGeneric_speciesMCMC_nfix(
                              name=name, params=params, M=x$M, n=n,
                              var_prior=var_prior,Nrep=Nrep, 
                              alpha=alpha, Rmax=Rmax, M_DM=NULL,
                              parallel = TRUE,M_max=M_max,seed0=x$seed)}  )
    filename = paste0(save_name_base,name,"_v_",var_prior,"_nfix_",trim_params,".Rdat")
    if(save_exp)
      save(ExpRes_list,file = filename)
  }
}

