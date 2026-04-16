# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[4] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/CriminalData")
setwd(wd)


# Functions ---------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
source("../../R/Rfunctions.R")
source("../../R/PFFAfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")

# Plot options ------------------------------------------------------------------
save_img = FALSE
width = 12; height = 6
cex.labels <- cex.lab <- 2
cex.axis <- 2
cex.legend <- 1.5
mycol = c("darkred","darkblue","lightblue","darkorange","darkgreen")
mycol2 = c("black","lightblue")
lgd_names = c("IBP","MBP","FB","Bounded","Unbounded")


# Load --------------------------------------------------------------------

load("RawDataInc.Rdat")
data = t(A)
n = ncol(A)
Kn = nrow(A)
N_j = rowSums(A)
names(N_j) = as.character(1:Kn)
Nj_ordered = sort(N_j, decreasing = TRUE)

seed = 34231
set.seed(seed)


# Options  --------------------------------------------------------
eps_grid = c(0.001, seq(0.1,0.3,length.out =  (34*5-1)) )
alpha = 0.05
M_max = 200
nstart = 10

seed0 = 4224
num_cores = 33 # <--- modify here
Nrep = 10 # <--- modify here

var_fct = 100
# Run) Mmax-based  --------------------------------------------------------
cat("\n Running stopping rule ... ")
res = SRinc_grid(eps_grid=eps_grid, data=data, nstart=nstart,
                 Nrep=Nrep, num_cores=num_cores, seed0=seed0,
                 alpha=alpha, var_fct=var_fct)
cat("done! Save and conclude \n")
save(res, file = "save/Mod3_Inc4People_SRMmax.Rdat")

# Run) Coverage-based  --------------------------------------------------------
# res_cov = SRabu_cov_grid( cov_grid, data, nstart, Nrep, num_cores, seed0)
# save(res_cov, file = "save/Mod2_Inc4Meetings_SRMcov.Rdat")

# Plot --------------------------------------------------------------------
stop_here = TRUE
ltype = 1
ygrid = eps_grid
if(!stop_here){
  load("save/Mod3_Inc4People_SRMmax.Rdat")
  # load("save/Mod2_Inc4Meetings_SRMcov.Rdat")
  
  res = lapply(res, function(x) {x[which(is.na(x))] = n; x} )
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
  
  if(save_img)
    pdf("img/Crim4features_StopR_grid.pdf", width = width, height = height)
  par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0), las = 1, cex = 2)
  plot(0,0,type = "n", main = "", ylab = "Nstop",
       xlim = range(eps_grid), ylim = c(0,n), 
       xlab = expression(epsilon) )
  # xlab = paste0(expression(epsilon)," / 1 - coverage") )
  for(i in 1:dim(res_all)[2]){
    points(y = res_all[2,i,], x = ygrid, 
           type = "l", lty = 1, 
           lwd = 3, col = mycol[i], pch = 16)
  }
  legend("bottomleft",lgd_names,
         fill = c(mycol),
         cex = cex.legend, bty = "n", border = NA)
  if(save_img)
    dev.off()
}










