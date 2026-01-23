# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/CriminalData/")
setwd(wd)

# Functions ---------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
source("../../R/Rfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# Colors ------------------------------------------------------------------


# Load --------------------------------------------------------------------

load("Locale.Rdat")
Locale_temp = Locale[-which(Locale == "OUT" | Locale == "MISS")]
data = Locale_temp
Nj_locale = table(Locale_temp)
Nj_locale = sort(Nj_locale, decreasing = TRUE)

n = sum(Nj_locale)
Kn = length(Nj_locale)


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
# res = SRabu_grid( eps_grid, data, nstart, Nrep, num_cores, seed0, alpha, M_max)
# save(res, file = "save/Mod1Abu_SRMmax.Rdat")

# Run) Coverage-based  --------------------------------------------------------
# res_cov = SRabu_cov_grid( cov_grid, data, nstart, Nrep, num_cores, seed0)
# save(res_cov, file = "save/Mod1Abu_SRcov.Rdat")

# Plot --------------------------------------------------------------------
stop_here = TRUE
ltype = c(1,1,1,2)
mycol = c("darkblue","darkred","darkorange","deeppink")
ygrids = vector("list",4)
ygrids[[1]]<-ygrids[[2]]<-ygrids[[3]]<-eps_grid
ygrids[[4]]<-(1-cov_grid)

if(!stop_here){
  load("save/Mod1Abu_SRMmax.Rdat")
  # load("save/Mod1Abu_SRcov.Rdat")
  
  res_list = lapply(res, function(x) apply(x,2,quantile,probs = c(0.025,0.5,0.975)))
  # res_cov_list = lapply(res_cov, function(x) apply(x,2,quantile,probs = c(0.025,0.5,0.975)))
  res_arr <- simplify2array(res_list) # 3 x 3 x length(eps_grid)
  # res_cov_arr <- simplify2array(res_cov_list) # 3 x 1 x length(eps_grid)
  # res_all <- array(
  #   do.call(cbind, lapply(seq_len(dim(res_arr)[3]), function(k)
  #     cbind(res_arr[,,k], res_cov_arr[,,k])
  #   )),
  #   dim = c(3, 4, 70)
  # )# 3 x 4 x length(eps_grid)
  res_all = res_arr
  
  par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0), las = 1)
  plot(0,0,type = "n", main = "", ylab = "Nstop",
       xlim = range(eps_grid), ylim = c(0,n), 
       # xlab = paste0(expression(epsilon)," / 1 - coverage") )
       xlab = paste0(expression(epsilon)) )
  for(i in 1:3){
    points(y = res_all[2,i,], x = ygrids[[i]], 
           type = "l", lty = ltype[i], 
           lwd = 3, col = mycol[i], pch = 16)
  }
  # abline(v = 0.05, lty = 2, col = "red")
  legend("bottomleft", c("FD","PYP","Freq"),#,"Cov"), 
         col = mycol[1:3], lty = ltype, lwd = 3)
  
}

