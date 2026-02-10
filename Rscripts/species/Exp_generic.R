# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[3] # <--- modify here
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
Exp_species_nfix_run = function(M,n,name,alpha = 0.05,Rmax = 100, M_max = 200, seed = 121321)
{
  source("../../R/Rfunctions.R")
  Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
  # From BinomialCIs
  source("../../../BinomialCIs/R/Rfunctions.R")
  Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")
  
  set.seed(seed)
  ptrue = sim_generic_species(name,M,n,alpha)
  # Define return object
  res_names = c("Mmax","Freq","PD","FDP")
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
  
  #e) Save results
  res[1,1] = Mmax
  res[1,2] = pain
  res[1,3] = ubpyp
  res[1,4] = ubFD
  
  # Return
  res
}

Exp_species_nfix = function(name,M,n,Nrep = 100, 
                            alpha = 0.05,Rmax = 100, 
                            M_max = 200, seed0 = 121321,
                            parallel = TRUE, num_cores = 5)
{
  cat("\n M = ",M,"\n")
  set.seed(seed0)
  seeds = sample(1:999999, size = Nrep)
  
  # Sequential case
  if(!parallel){
    res_list = lapply(seeds, function(seed) Exp_species_nfix_run(M=M, n=n, alpha=alpha, 
                                                                 Rmax=Rmax, M_max=M_max, 
                                                                 seed=seed, name=name) )
  }else{
    ## Parallel case
    cluster <- makeCluster(num_cores, type = "SOCK")
    doSNOW::registerDoSNOW(cluster)
    clusterExport(cluster, list("alpha"), envir = environment())
    res_list = parLapply( cl = cluster, x = seeds,
                          fun = Exp_species_nfix_run,
                          M=M, n=n, alpha=alpha, name=name,
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
mycol = c("darkorange","darkred","darkblue","black")
mycol_band = c("lightsalmon","tomato","lightblue","grey85")
lgd_names = c("Freq","PD","FDP","Oracle")

# Options -----------------------------------------------------------------
experiments = list("WorstUnif")
alpha <- alfa <- 0.05
Nrep = 50 
n = 500
Rmax = 100; RmaxFD = 50
Mmin_grid = 50; Mmax_grid = 1000
Mgrid = seq(Mmin_grid,Mmax_grid,by = 50); LMgrid = length(Mgrid)
Nexp = length(Mgrid)
M_max = 200

seed0 = 42

save_exp = TRUE # <---
save_name_base = paste0("save/Species_") 
img_fld = paste0("img/") 
# n fix -----------------------------------------------------------------

ii = 1
name = experiments[[ii]]

num_cores = 35 # <---
ExpRes_list = lapply(Mgrid, function(M) Exp_species_nfix(name=name, M=M, n=n,
                                                         Nrep=Nrep, alpha=alpha, Rmax=Rmax, 
                                                         M_max=M_max,seed0=seed0) )


## Oracle 
ma = find_ma_worstunif(n,alpha)
oracle = oracle_worst_uniform(n = n, ma = ma, alpha = alpha)

ExpRes_list = lapply(ExpRes_list, function(x){ y = matrix(oracle, nrow = Nrep, ncol = 1); colnames(y) = "oracle"; cbind(x,y) } )

filename = paste0(save_name_base,name,"_nfix",".Rdat")
if(save_exp)
  save(ExpRes_list, file = filename)


# Analysis ----------------------------------------------------------------

stop_here = TRUE

if(!stop_here){
  # Coverage ----------------------------------------------------------------
  cov_list = lapply(ExpRes_list, function(x){
    cov_mat = matrix(NA,nrow = 1, ncol = 4)
    cov_mat[1,1] = length(which(x[,1] <= x[,2]))
    cov_mat[1,2] = length(which(x[,1] <= x[,3]))
    cov_mat[1,3] = length(which(x[,1] <= x[,4]))
    cov_mat[1,4] = length(which(x[,1] <= x[,5]))
    cov_mat/nrow(x)
  })
  cov_mat = do.call(rbind,cov_list)
  cov_mat
  
  
  # Plot --------------------------------------------------------------------
  ExpRes_qnt = lapply(ExpRes_list, function(x) apply(x,2,quantile,probs = c(0.025,0.5,0.975)) )
  ExpRes_qnt <- simplify2array(ExpRes_qnt) # 3 x 5 x LMgrid
  
  ## axis labels
  ymax = (11/10) * max(ExpRes_qnt); ymin = 0
  ymax = 15*1e-3; ymin = 11.5*1e-3
  ylim_plot = c(ymin,ymax)
  ypos = seq(ymin,ymax,length.out = 5)
  ylabs = as.character(round(ypos*1e3,0))
  xmax = max(Mgrid); xmin = min(Mgrid)
  xlim_plot = c(xmin,xmax)
  xpos = Mgrid
  xlabs = as.character(Mgrid)
  
  
  img_name = paste0(img_fld,name,"_nfix",".pdf")
  if(save_img)
    pdf(img_name, width = width, height = height)
  par( mfrow = c(1,1), mar = c(3.5,4,1,1), mgp=c(2.5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
  plot(0,0,  yaxt = "n", xaxt = "n",
       xlab = "", ylab = "1000 * bound",
       xlim = xlim_plot , ylim = ylim_plot, 
       main = paste0(" "),
       type = "n")
  grid(lty = 1,lwd = 1, col = "gray90" )
  axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
  axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
  mtext("M", side = 1, line = 2.5, cex = cex.axis)
  for(ii in 2:dim(ExpRes_qnt)[2]){
    points( x = Mgrid, y = ExpRes_qnt[2,ii,], 
            type = "l", lwd = 5, col = mycol[ii-1] ) 
    # polygon( c(Mgrid, rev(Mgrid)),
    #          c(ExpRes_qnt[1,ii,], rev(ExpRes_qnt[3,ii,])),
    #          col = mycol_band[ii-1],
    #          border = NA) # plot in-sample bands
  }
  legend("right",lgd_names,
         fill = mycol, cex = cex.legend, bty = "n", border = NA)
  if(save_img)
    dev.off()
  
  
  # Brutta ------------------------------------------------------------------
  
  
  # #a) Oracle
  # if(name == "WorstUnif"){
  #   ma = find_ma_worstunif(n,alpha)
  #   oracle = oracle_worst_uniform(n = n, ma = ma, alpha = alpha)
  # }else{
  #   # Oracle
  #   Mmax = rep(NA,Bor)
  #   for(b in 1:Bor){
  #     # Gen. data
  #     
  #     if(Kn == M){
  #       Mmax[b] = 0
  #     }else{
  #       idx_unobs = which(n_i == 0)
  #       Mmax[b] = max(ptrue[idx_unobs])
  #     }
  #   }
  #   oracle = quantile(Mmax, probs = 1-alfa)
  # }
  
}
