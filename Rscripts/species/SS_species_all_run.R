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
params_zipfs  = list(0.9,1.02,2,3)
params_geom   = list(0.85,0.9,0.95)
params_unif   = list(NA)
params_negbin = list(c(1,0.003),c(5,0.003),c(1,0.01))
experiments   = list("Zipfs"   = params_zipfs,
                     "Geom"    = params_geom,
                     "Uniform" = params_unif,
                     "NegBin"  = params_negbin)

alpha <- alfa <- 0.05
num_cores = 33 # <---
Nrep = 500 # <---
Rmax = 100; RmaxFD = 50
M_max = 200
seed0 = 42
set.seed(seed0)
var_prior = 10
parallel = TRUE

# n fix -----------------------------------------------------------------
# n = 500
# Mmin_grid = 50; Mmax_grid = 1000
# Mgrid = seq(Mmin_grid,Mmax_grid,by = 50); LMgrid = length(Mgrid)
# Nexp = length(Mgrid)
# seeds = sample(1:999999, size = Nexp)
# 
# save_exp = TRUE # <---
# save_name_base = paste0("save/SS_species_nfixMAP_")
# img_fld = paste0("img/")
# igrid = c(1:4)
# ii = 4
# for(ii in igrid){
# 
#   name = names(experiments)[ii]
#   Ncases = length(experiments[[ii]])
#   jj = 1
#   for(jj in 1:Ncases){
#     cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
#     params = experiments[[ii]][[jj]]
#     trim_params = sapply(params, get_first3digits, 4)
#     if(length(trim_params)>1)
#       trim_params = paste0(trim_params[1],"_",trim_params[2])
# 
#     run_obj <- Map(function(m, s) list(M = m, seed = s),Mgrid, seeds)
# 
#     ExpRes_list = lapply( run_obj,
#                           function(x){
#                             SS_species(
#                               name=name, params=params, M=x$M, n=n,
#                               var_prior=var_prior,Nrep=Nrep,
#                               alpha=alpha, Rmax=Rmax, M_DM=NULL,
#                               parallel = parallel,M_max=M_max,seed0=x$seed)}  )
#     filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
#     if(save_exp)
#       save(ExpRes_list,file = filename)
#   }
# }
# 

# M fix -----------------------------------------------------------------
M = 1000
Nmin_grid = 100; Nmax_grid = 5000
Ngrid = seq(Nmin_grid,Nmax_grid,by = 500); LNgrid = length(Ngrid)
Nexp = length(Ngrid)
seeds = sample(1:999999, size = Nexp)

save_exp = TRUE # <---
save_name_base = paste0("save/SS_species_MfixMAP_") 
img_fld = paste0("img/") 
igrid = c(1:2)
ii = 1
for(ii in igrid){
  
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 4
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = sapply(params, get_first3digits, 4)
    if(length(trim_params)>1)
      trim_params = paste0(trim_params[1],"_",trim_params[2])
    
    run_obj <- Map(function(n, s) list(N = n, seed = s),Ngrid, seeds)
    
    ExpRes_list = lapply( run_obj,
                          function(x){
                            SS_species(
                              name=name, params=params, M=M, n=x$N,
                              var_prior=var_prior,Nrep=Nrep, 
                              alpha=alpha, Rmax=Rmax, M_DM=NULL,
                              parallel = parallel,M_max=M_max,seed0=x$seed)}  )
    filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
    if(save_exp)
      save(ExpRes_list,file = filename)
  }
}
