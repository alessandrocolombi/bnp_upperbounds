# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/species")
setwd(wd)

# Functions ---------------------------------------------------------------
source("../../R/Rfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")

# Plot options ------------------------------------------------------------------
save_img = FALSE
width = 12; height = 6
cex.labels <- cex.lab <- 2
cex.axis <- 2
cex.legend <- 1.5
# Custom function ---------------------------------------------------------

# Example -----------------------------------------------------------------
alpha = 0.05
M = 500
n = 1000

ptrue <- worst_uniform(M = M, n = n, alpha = alpha)
ptrue


# Ami's paper example -----------------------------------------------------
alpha <- alfa <- 0.05
n = 1000
Rmax = 100
seed = 23131
set.seed(seed)

Mmin = 10; Mmax = 1000; LMgrid = 20
Mgrid = round(seq(Mmin,Mmax,length.out = LMgrid))
Nexp = LMgrid

pain = ub_pain(n = n, Rmax = Rmax, alfa = alfa)
pain = rep( pain, Nexp )
ma = find_ma_worstunif(n,alpha)
oracle = oracle_worst_uniform(n = n, ma = ma, alpha = alpha)
oracle = rep(oracle, LMgrid)
ii = 15
for(ii in seq_along(Mgrid)){
  M = Mgrid[ii]
  ptrue <- worst_uniform(M = M, n = n, alpha = alpha)
  
  # Oracle
  # Mmax = rep(NA,Bor)
  # b = 1
  # for(b in 1:Bor){
  #   # Gen. data
  #   data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
  #   n_i = tabulate(data, nbins = M)
  #   idx_obs = which(n_i > 0)
  #   Kn = length(idx_obs)
  #   data_obs = n_i[idx_obs]
  #   if(Kn == M){
  #     Mmax[b] = 0
  #   }else{
  #     idx_unobs = which(n_i == 0)
  #     Mmax[b] = max(ptrue[idx_unobs])
  #   }
  # }
  # oracle[ii] = quantile(Mmax, probs = 1-alfa)
}



# Plot --------------------------------------------------------------------

ymax = max(oracle,pain);ymin = 0
ylabs = round(seq(ymin*1e3,ymax*1e3,length.out = 5),1)

img_name = "img/WorstCase_nfix"

if(save_img)
  pdf(img_name, width = width, height = height)
par( mfrow = c(1,1), mar = c(3.5,5,1,1), mgp=c(3.5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,  yaxt = "n", xaxt = "n",
     xlab = "", ylab = "1000 * bound",
     xlim = c(min(Mgrid),max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = cex.axis )
axis(1, cex.axis = cex.axis);mtext("M", side = 1, line = 2.5, cex = cex.axis)
points( x = Mgrid, y = oracle, 
        type = "l", 
        lwd = 5, pch = 16, lty = 1,
        col = "black" ) 
points( x = Mgrid, y = pain,
        type = "l",
        lwd = 5, pch = 16, lty = 1,
        col = "darkorange" )
legend("bottomright", legend = c("Oracle","Freq."), 
       fill = c("black","darkorange"), 
       cex = cex.legend, bty = "n", border = NA)
if(save_img)
  dev.off()







