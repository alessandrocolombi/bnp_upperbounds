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
Nrep = 100
n = 200
M = 2000 
M_max = 200
seed = 232131
set.seed(seed)
Mmax = M;
AccCrv_length = 50
ii = 1; name = names(experiments)[[ii]]
jj = 1
params = experiments[[ii]][[jj]]
trim_params = sapply(params, get_first3digits, 4)
if(length(trim_params)>1)
  trim_params = paste0(trim_params[1],"_",trim_params[2])

ptrue = sim_generic_species(name,M,params)
ptrue = sort(ptrue, decreasing = TRUE)
ptrue_mat = matrix(nrow = Nrep,ncol = M)
ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)
data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
Nj_mat = apply(data_mat, 2, tabulate, nbins = M) # M x Nrep

# Acc curve ---------------------------------------------------------------

n_reorderings = 10
orderings = lapply(1:n_reorderings, function(j) sample(1:n, size = n))
AccCurves_list = lapply(orderings, function(o) data_mat[o,1])
AccCurves_list = lapply(AccCurves_list, function(x){
  sapply(1:length(x), function(j) { y = x[1:j]; length(table(y)) })
})
AccCurves_mat <- do.call(rbind, AccCurves_list) #(n_reorderings x n_people)
AccCurves_res = apply(AccCurves_mat, 2, quantile, probs = c(0.025,0.5,0.975))

# Plot
ngrid = 1:n

par(mfrow = c(1,1), mar = c(2.5,4.5,0,1), mgp=c(3,0.5,0), bty = "l", cex = 1)
plot(x = 0, y = 0, type = "n",
     main = " ", xlab = "n", ylab = "#symbols",
     ylim = c(0,max(AccCurves_res)+1),
     xlim = c(0,n+1),
     pch = 1) # init plot
polygon( c(ngrid, rev(ngrid)),
         c(AccCurves_res[1,], rev(AccCurves_res[3,])),
         col = "grey75",
         border = NA) # plot in-sample bands
points(x = ngrid, y = AccCurves_res[2,], type = "l", lwd = 3) # plot mean obs


# Fit -----------------------------------------------------
var_prior = 10000
Rmax = 100
alpha <- alfa <- 0.05
seed = 213271313
fit = t(apply(data_mat, 2, function(data){
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
  #b) Run analysis
  UB_fit(n=n,Kn=Kn,n_i=n_i,data_obs=data_obs,
         Mmax=Mmax,M=M,var_prior=var_prior,Rmax=Rmax,
         alpha=alpha,useMAP = TRUE,M_max=M_max,seed=seed)
}))
colnames(fit) = c("Mmax","Freq","PD","FDP","DirMulti", "sigma","theta","gamma_FDP","Lambda_FDP","gamma_DM", "Kn")

fit[,c(1,3,4,8,9,11)]

# Mstar -------------------------------------------------------------------


ExpMstar = apply(fit, 1, function(x){
  gamma = x[8]; Lambda = x[9]; Kn = x[11]
  Mstar_ub = 5000
  ## Mstar
  logExpMstar = compute_logV( Kn+1, n, gamma, Lambda, M_max ) - 
    compute_logV( Kn,   n, gamma, Lambda, M_max )
  exp(logExpMstar) # Expected value
})
fit = cbind(fit, ExpMstar)


fit[,c(9,11,12)]


# Boxplots for Lambda, Kn, ExpMstar --------------------------------------

par( mfrow = c(1,3), mar = c(1,4,1,1), mgp=c(2.75,0.5,0), bty = "l", 
     las = 1, cex = 1 )
for(ij in c(9,11,12)) {
  boxplot(fit[,ij],
          ylab = colnames(fit)[ij], xlab = "",
          col = "grey85", border = "grey30")
}


# Coverage ----------------------------------------------------------------

x = fit[,1:5]
counts <- colSums(x[, -1, drop = FALSE] >= x[, 1])
counts / Nrep


apply(fit[,1:5], 2, quantile, probs = c(0.025,0.5,0.975))
