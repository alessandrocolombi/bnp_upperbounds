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

# Options -----------------------------------------------------------------
experiments = list("WorstUnif")
alpha <- alfa <- 0.05
Nrep = 500 
n = 500
Rmax = 100; RmaxFD = 50
Mmin_grid = 50; Mmax_grid = 500
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
ii = 1
name = experiments[[ii]]
ma = find_ma_worstunif(n,alpha)
q = 1-p_all_seen_uniform(n,ma)
trim_q = get_first3digits(q,4)

# Analysis ----------------------------------------------------------------
filename = paste0(save_name_base,name,"_n",n,"_q",trim_q,".Rdat")
load(filename)

  
# Coverage ----------------------------------------------------------------
cov_mat <- do.call(rbind, lapply(ExpRes_list, function(x) {
  n <- nrow(x)
  if (n == 0) return(rep(NA_real_, ncol(x) - 1))
  counts <- colSums(x[, -1, drop = FALSE] >= x[, 1])
  counts / n
}))
colnames(cov_mat) <- colnames(ExpRes_list[[1]])[-1]
cov_mat

# Plot - mediana --------------------------------------------------------------------
ExpRes_qnt = lapply(ExpRes_list, function(x) apply(x[,c(2:6)],2,quantile,probs = c(0.025,0.5,0.975)) )
ExpRes_qnt <- simplify2array(ExpRes_qnt) # 3 x 5 x LMgrid
  
## axis labels
ymax = (11/10) * max(ExpRes_qnt); ymin = 10.5*1e-3
# ymax = 15*1e-3; ymin = 8.5*1e-3
ylim_plot = c(ymin,ymax)
ypos = seq(ymin,ymax,length.out = 5)
ylabs = as.character(round(ypos*1e3,0))
xmax = max(Mgrid); xmin = min(Mgrid)
xlim_plot = c(0,xmax)
xpos = Mgrid
xlabs = as.character(Mgrid)
  
  
img_name = paste0(img_fld,name,"_n",n,"_q",trim_q,".pdf")
if(save_img)
  pdf(img_name, width = width, height = height)
par( mfrow = c(1,1), mar = c(3.5,4.25,1,1), mgp=c(2.75,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,  yaxt = "n", xaxt = "n",
       xlab = "", ylab = "1000 * bound",
       xlim = xlim_plot , ylim = ylim_plot, 
       main = paste0(" "),
       type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
mtext("M", side = 1, line = 2.5, cex = cex.axis)
abline(h = 1/ma, lty = 1, lwd = 2, col = "black")
for(ii in 1:dim(ExpRes_qnt)[2]){
  points( x = Mgrid, y = ExpRes_qnt[2,ii,], 
          type = "l", lwd = 5, col = mycol[ii] )
  points( x = Mgrid, y = ExpRes_qnt[1,ii,], 
          type = "l", lwd = 3, lty = 2, col = mycol[ii] )
  points( x = Mgrid, y = ExpRes_qnt[3,ii,], 
          type = "l", lwd = 3, lty = 2, col = mycol[ii] )
}
legend("topleft",lgd_names,
         fill = mycol, 
       cex = cex.legend, bty = "n", border = NA)
if(save_img)
  dev.off()
  

# Dir-Multi-m analysis -------------------------------------------------------------


DirMulti_Ma_list = lapply(ExpRes_list, function(x) {
    not_all_seen = which(x[,1] > 0)
    y = x[not_all_seen,5]
    y
    # quantile(y, probs = c(0.025,0.5,0.975))
    # quantile(y, probs = c(0.025,0.5,0.975))
})
DirMulti_Ma_qnt = lapply(ExpRes_list, function(x) {
  not_all_seen = which(x[,1] > 0)
  y = x[not_all_seen,5]
  quantile(y, probs = c(0.025,0.5,0.975))
})
DirMulti_Ma_qnt = do.call(cbind,DirMulti_Ma_qnt)

min(sapply(DirMulti_Ma_list, min))
 
  ## axis labels
  ymax = (11/10) * max(sapply(DirMulti_Ma_list,max)); ymin = 0
  ylim_plot = c(ymin,ymax)
  ypos = seq(ymin,ymax,length.out = 5)
  ylabs = as.character(round(ypos*1e3,0))
  xmax = max(Mgrid); xmin = min(Mgrid)
  xlim_plot = c(xmin,xmax)
  xpos = Mgrid
  xlabs = as.character(Mgrid)
  
  
  img_name = paste0(img_fld,name,"_nfix_DirMultiMa",".pdf")
  if(save_img)
    pdf(img_name, width = width, height = height)
  par( mfrow = c(1,1), mar = c(3.5,4.25,1,1), mgp=c(2.75,1,0), bty = "l", las = 1, cex.lab = cex.lab )
  plot(0,0,  yaxt = "n", xaxt = "n",
       xlab = "", ylab = "1000 * bound",
       xlim = xlim_plot , ylim = ylim_plot, 
       main = paste0(" "),
       type = "n")
  grid(lty = 1,lwd = 1, col = "gray90" )
  axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
  axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
  mtext("M", side = 1, line = 2.5, cex = cex.axis)
  abline(v = ma, lty = 2, lwd = 2, col = "red")
  abline(h = 1/ma, lty = 3, lwd = 2, col = "black")
  for(ii in 1:length(DirMulti_Ma_list)){
    boxplot(DirMulti_Ma_list[[ii]], add = TRUE,
            at = Mgrid[ii], yaxt = "n")
  }

# Brutta ------------------------------------------------------------------

  
prova = lapply(ExpRes_list[2:10], function(x) {
  not_all_seen = which(x[,1] > 0)
  y = x[not_all_seen,c(3:6)]
})
prova  

sopra = lapply(prova, function(x){
  length(which(x[,2]< (1/ma)))/nrow(x)
})
sopra
