# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/species/")
setwd(wd)

# Functions ---------------------------------------------------------------
source("../../R/Rfunctions.R")
source("../../R/pyp_utils_fast.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")

# Funzioni ----------------------------------------------------------------

# Common params -----------------------------------------------------------
seed = 42
set.seed(seed)
alfa <- alpha <- 0.05
n = 1000

# Sim. data from model -------------------------------------------------------
BB = 100

sigma_true = 0.25; theta_true = 50
xmax = 500

sim_PD_mat = r_SB(BB,xmax, sigma_true, theta_true,seed) # BB x xmax
sim_PD_mat = apply(sim_PD_mat, 1, sort, decreasing = TRUE) # xmax x BB
colSums(sim_PD_mat)
qnt_PD = apply(sim_PD_mat, 1, quantile, probs = c(0.025,0.5,0.975)) # 3 x xmax

## Plot data
ymax = max(qnt_PD) * 1.02
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,100), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_PD[1,1:xmax], rev(qnt_PD[3,1:xmax])),
         col = "grey85",
         border = NA) # plot in-sample bands

# Gen. all datasets
counts <- vapply( seq_len(BB), function(b) as.integer(rmultinom(1, size = n, prob = sim_PD_mat[, b])),
                  integer(xmax) ) # xmax x BB
dim(counts)
colSums(counts) 

## Fit and check Ale -----------------------------------------------------------
# All datasets
start_params <- c(alpha = 0.5, theta = 1)
res_vec = apply(counts, 2, function(Nj){
  Nj = Nj[which(Nj > 0)]
  Kn = length(Nj)
  fit <- tryCatch(optim(par = start_params, fn = llik_pyp, 
                        n = n, Kn = Kn, data_obs = Nj, # extra parameters
                        method = "L-BFGS-B",
                        lower = c(0, -1), upper = c(1-1e-10, Inf)) ,
                  error=function(e) NA)
  list("sigma_mle" = fit$par[1],
       "theta_mle" = fit$par[2],
       "fval" = fit$value)
})

res_mat = matrix(0,nrow = BB, ncol = 0)
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$sigma_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$theta_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$fval) )
colnames(res_mat) <- names(res_vec[[1]])

round(colMeans(res_mat),2)

fitcheck_all = apply(res_mat, 1, function(x) {
  r_SB(1,xmax, x[1], x[2], seed)
})
fitcheck_all = apply(fitcheck_all, 2, sort, decreasing = TRUE) 
fitcheck_qnt = apply(fitcheck_all, 1, quantile, probs = c(0.025,0.5,0.975))

## Plot data
ymax = max(max(qnt_PD),max(fitcheck_qnt)) * 1.02
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,100), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1)
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_PD[1,1:xmax], rev(qnt_PD[3,1:xmax])),
         col = "grey85",
         border = NA) # plot in-sample bands
polygon( c(1:xmax, rev(1:xmax)),
         c(fitcheck_qnt[1,1:xmax], rev(fitcheck_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands


## Fit and check Mario -----------------------------------------------------------


# Fit PYP parameters for each dataset (each column of counts)
fit_list <- apply(counts, 2, function(Nj) {
  sizes <- Nj[Nj > 0]  # observed species counts
  fit_pyp_profile(
    sizes,
    refine = TRUE,
    verbose = FALSE
  )
})

# Convert to matrix: sigma(alpha), theta, loglik
res_mat <- cbind(
  sigma_mle = sapply(fit_list, `[[`, "alpha"),
  theta_mle = sapply(fit_list, `[[`, "theta"),
  loglik    = sapply(fit_list, `[[`, "loglik"),
  converged = sapply(fit_list, `[[`, "converged")
)

round(colMeans(res_mat[, c("sigma_mle","theta_mle","loglik")]), 2)
table(res_mat[, "converged"])



fitcheck_all = apply(res_mat, 1, function(x) {
  r_SB(1,xmax, x[1], x[2], seed)
})
fitcheck_all = apply(fitcheck_all, 2, sort, decreasing = TRUE) 
fitcheck_qnt = apply(fitcheck_all, 1, quantile, probs = c(0.025,0.5,0.975))

## Plot data
ymax = max(max(qnt_PD),max(fitcheck_qnt)) * 1.02
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,100), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1)
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_PD[1,1:xmax], rev(qnt_PD[3,1:xmax])),
         col = "grey85",
         border = NA) # plot in-sample bands
polygon( c(1:xmax, rev(1:xmax)),
         c(fitcheck_qnt[1,1:xmax], rev(fitcheck_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands






# Worst Uniform -----------------------------------------------------------
M = 200
alpha = 0.05
Rmax = 100
n = 500 
name = "WorstUnif"
M_max = 200
Nrep = 100

params = c(n,alpha)
ptrue = sim_generic_species(name,M,params)
ptrue = sort(ptrue, decreasing = TRUE)
ptrue_mat = matrix(nrow = Nrep,ncol = M)
ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)
data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
Nj_mat = apply(data_mat, 2, tabulate, nbins = M) # M x Nrep


## Stima Ale
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
colnames(EB_WorstUnif) = c("sigma","theta","LambdaFDP","gammaFDP","gammaDM","Kn",
                           "Mmax","Freq",
                           "UB_PD","UB_FDP","UB_DM","pmin_obs")


EB_WorstUnif[,1];EB_WorstUnif[,2];


## Stima Mario
EB_WorstUnif2 = t(apply(data_mat, 2, function(data){
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
  fit = fit_pyp_profile(Kn,refine = TRUE,verbose = FALSE)
  res[1,1] = fit$alpha
  res[1,2] = fit$theta
  
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
colnames(EB_WorstUnif2) = c("sigma","theta","LambdaFDP","gammaFDP","gammaDM","Kn",
                           "Mmax","Freq",
                           "UB_PD","UB_FDP","UB_DM","pmin_obs")


EB_WorstUnif2[,1];EB_WorstUnif2[,2]
EB_WorstUnif2[,9]
