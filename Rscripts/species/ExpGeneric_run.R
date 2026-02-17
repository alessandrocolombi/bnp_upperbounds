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
ExpGeneric_species_nfix_run = function(M, n, name, params,
                                       alpha = 0.05,Rmax = 100, M_max = 200, M_DM = NULL, seed = 121321)
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
  res_names = c("Mmax","Freq","PD","FDP","DirMulti")
  res = matrix(NA,nrow = 1, ncol = length(res_names))
  colnames(res) = res_names
  #a) Generate data
  data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
  n_i = tabulate(data, nbins = M)
  idx_obs = which(n_i > 0)
  Kn = length(idx_obs)
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
  #c) PD
  # Param. estimation (PYP)
  start_params <- c(alpha = 0.5, theta = 1)
  fit <- optim(par = start_params, fn = llik_pyp, 
               n = n, Kn = Kn, data_obs = data_obs, # extra parameters
               method = "L-BFGS-B",
               lower = c(0, -1), upper = c(1-1e-10, Inf)) 
  alpha_mle = fit$par[1]
  theta_mle = fit$par[2]
  
  # Upper bound (PYP)
  ubpyp = exp(compute_log_UBMarkov( Rmax, alpha_mle, theta_mle, Kn, n, alpha ))
  ubpyp = min(ubpyp,1)
  #d) FDP
  # Param. estimation (FD)
  start_params <- c(gamma = 0.1, Lambda = Kn)
  fit <- optim(par = start_params, fn = llik_FD, 
               n = n, Kn = Kn, data_obs = data_obs, M_max = M_max,# extra parameters
               method = "L-BFGS-B",
               lower = c(1e-5, 1e-5), upper = c(Inf, Inf)) 
  gamma_mle = fit$par[1]
  Lambda_mle = fit$par[2]
  # Upper bound (FD)
  ubFD = exp(compute_log_UBMarkov_FD( Rmax, gamma_mle, Lambda_mle, Kn, n, alpha, M_max ))
  ubFD = min(ubFD,1)
  #f) Dirichlet-Multinomial --> M = M_DM
  M_DM = ifelse(is.null(M_DM), M, M_DM)
  # Param. estimation (DirMulti)
  start_params <- c(gamma = 0.5)
  fit <- optim(par = start_params, fn = llik_DirMult,
               n = n, M = M_DM, data = n_i, # extra parameters
               method = "L-BFGS-B",
               lower = c(1e-10), upper = c(Inf))
  gamma_mle = fit$par[1]
  # Upper bound (DirMulti)
  ubDM = exp(compute_log_UB_DirMulti( Rmax, gamma_mle, M, Kn, n, alpha ))
  ubDM = min(ubDM,1)
  #f) Save results
  res[1,1] = Mmax
  res[1,2] = pain
  res[1,3] = ubpyp
  res[1,4] = ubFD
  res[1,5] = ubDM
  # Return
  res
}

ExpGeneric_species_nfix = function(name,M,n,params,Nrep = 100, 
                                   alpha = 0.05,Rmax = 100, 
                                   M_max = 200, seed0 = 121321,
                                   M_DM = NULL,
                                   parallel = TRUE, num_cores = 5)
{
  cat("\n M = ",M,"\n")
  set.seed(seed0)
  seeds = sample(1:999999, size = Nrep)
  
  # Sequential case
  if(!parallel){
    res_list = lapply(seeds, function(seed) ExpGeneric_species_nfix_run(
                                              M=M, n=n, params=params, 
                                              alpha=alpha, 
                                              Rmax=Rmax, M_max=M_max,
                                              M_DM=M_DM,
                                              seed=seed, name=name) )
  }else{
    ## Parallel case
    cluster <- makeCluster(num_cores, type = "SOCK")
    doSNOW::registerDoSNOW(cluster)
    clusterExport(cluster, list("alpha"), envir = environment())
    res_list = parLapply( cl = cluster, x = seeds,
                          fun = ExpGeneric_species_nfix_run,
                          M=M, n=n, alpha=alpha, name=name, params=params, 
                          M_DM=M_DM,
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
num_cores = 25 # <---
Nrep = 5000 # <---
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
save_name_base = paste0("save/Species_") 
img_fld = paste0("img/") 
# n fix -----------------------------------------------------------------

igrid = c(4)
ii = 1
for(ii in igrid){
  
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 3
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = sapply(params, get_first3digits, 4)
    if(length(trim_params)>1)
      trim_params = paste0(trim_params[1],"_",trim_params[2])
    
    run_obj <- Map(function(m, s) list(M = m, seed = s),Mgrid, seeds)
    
    ExpRes_list = lapply( run_obj,
                          function(x){
                            ExpGeneric_species_nfix(
                              name=name, params=params, M=x$M, n=n,
                              Nrep=Nrep, alpha=alpha, Rmax=Rmax,
                              M_DM=NULL,
                              parallel = TRUE,
                              M_max=M_max,seed0=x$seed)}  )
    filename = paste0(save_name_base,name,"_nfix_",trim_params,".Rdat")
    if(save_exp)
      save(ExpRes_list,file = filename)
  }
}

