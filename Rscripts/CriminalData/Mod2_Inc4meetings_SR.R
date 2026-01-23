# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[3] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/CriminalData/")
setwd(wd)

# Functions ---------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
source("../../R/Rfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")

# Colors ------------------------------------------------------------------


# Load --------------------------------------------------------------------

load("RawDataInc.Rdat")
data = A
n = nrow(A)
Kn = ncol(A)
N_j = colSums(A)
names(N_j) = as.character(1:Kn)

seed = 34231
set.seed(seed)


# Options  --------------------------------------------------------
eps_grid = seq(0.001,0.2,length.out =  34*5)
cov_grid = 1 - eps_grid
alpha = 0.05
M_max = 200
nstart = 10

seed0 = 4224
num_cores = 34
Nrep = 50

# Run) Mmax-based  --------------------------------------------------------
res = SRinc_grid( eps_grid, data, nstart, Nrep, num_cores, seed0, alpha)
save(res, file = "save/Mod2_Inc4Meetings_SRMmax.Rdat")

# Run) Coverage-based  --------------------------------------------------------
# res_cov = SRabu_cov_grid( cov_grid, data, nstart, Nrep, num_cores, seed0)
# save(res_cov, file = "save/Mod2_Inc4Meetings_SRMcov.Rdat")

# Plot --------------------------------------------------------------------
stop_here = TRUE
ltype = c(1,1,1,2)
mycol = c("darkblue","darkred","darkorange","deeppink")
ygrids = vector("list",4)
ygrids[[1]]<-ygrids[[2]]<-ygrids[[3]]<-eps_grid
ygrids[[4]]<-(1-cov_grid)

if(!stop_here){
  load("save/Mod2_Inc4Meetings_SRMmax.Rdat")
  load("save/Mod2_Inc4Meetings_SRMcov.Rdat")
  
  res_list = lapply(res, function(x) apply(x,2,quantile,probs = c(0.025,0.5,0.975)))
  res_cov_list = lapply(res_cov, function(x) apply(x,2,quantile,probs = c(0.025,0.5,0.975)))
  res_arr <- simplify2array(res_list) # 3 x 3 x length(eps_grid)
  res_cov_arr <- simplify2array(res_cov_list) # 3 x 1 x length(eps_grid)
  res_all <- array(
    do.call(cbind, lapply(seq_len(dim(res_arr)[3]), function(k)
      cbind(res_arr[,,k], res_cov_arr[,,k])
    )),
    dim = c(3, 4, 70)
  )# 3 x 4 x length(eps_grid)
  
  par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
  plot(0,0,type = "n", main = "", ylab = "Nstop",
       xlim = range(eps_grid), ylim = c(0,n), 
       xlab = paste0(expression(epsilon)," / 1 - coverage") )
  for(i in 1:4){
    points(y = res_all[2,i,], x = ygrids[[i]], 
           type = "l", lty = ltype[i], 
           lwd = 3, col = mycol[i], pch = 16)
  }
  # abline(v = 0.05, lty = 2, col = "red")
  legend("bottomleft", c("FD","PYP","Freq","Cov"), 
         col = mycol, lty = ltype, lwd = 3)
  
}

