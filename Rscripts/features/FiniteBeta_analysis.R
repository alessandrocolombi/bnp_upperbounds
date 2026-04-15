
# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/features")
setwd(wd)

# Functions ---------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))

source("../../R/Rfunctions.R")
source("../../R/PFFAfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")



# Plot options ------------------------------------------------------------------
save_img = TRUE
width = 12; height = 6
cex.labels <- cex.lab <- 2
cex.axis <- 2
cex.legend <- 1.5
mycol = c("darkgreen","darkorange","darkred","darkblue","lightblue")

lgd_names = c("oracle","Bounded","Unbounded","IBP","MBP","FB")


# Ub shape ----------------------------------------------------------------

Rmax = 100
n = 500
M = 100
Kn = 90
alpha <- alfa <- 0.05
tol = sqrt(.Machine$double.eps) 

mean_grid = seq(1e-5,1-1e-5,length.out = 1000)
kappa_grid = c(1+2*tol,2,10,50,100)
UB_mat = matrix(-1, nrow = length(mean_grid), ncol = length(kappa_grid))

for(i in seq_along(mean_grid)){
  vmax = mean_grid[i] * (1 - mean_grid[i])
  for(j in seq_along(kappa_grid)){
    v = vmax/kappa_grid[j] 
    ab = compute_ab_beta(mean_grid[i],v)
    UB_mat[i,j] = min(1, exp(compute_log_UBMarkov_FB(Rmax,ab[1],ab[2],n,Kn,M,alpha) ) )
  }
}

mycol_ub = hcl.colors(n = length(kappa_grid))
## axis labels
ymax = (11/10) * max(UB_mat); ymin = (10/11) * min(UB_mat)
ylim_plot = c(ymin,ymax)
ypos = seq(ymin,ymax,length.out = 5)
ylabs = as.character(round(ypos,2))
xmax = max(mean_grid); xmin = min(mean_grid)
xpos = seq(xmin,xmax,length.out = 10)
xlabs = as.character(floor(xpos*100))
xlim_plot = c(xmin,xmax)

par( mfrow = c(1,1), mar = c(3.5,6,1,1), mgp=c(4.5,1,0), bty = "l", las = 1, cex.lab = 2 )
plot(0,0,  yaxt = "n", xaxt = "n",
     xlab = "", ylab = "bound",
     xlim = xlim_plot , ylim = ylim_plot, 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
mtext("mean * 100", side = 1, line = 2.5, cex = cex.axis)
for(j in seq_along(kappa_grid)){
  points( x = mean_grid, y = UB_mat[,j], 
          type = "l", lwd = 5, col = mycol_ub[j] )
}




# zipfs fail --------------------------------------------------------------

M = 5000
n = 500
Kn = 4998
oracle = 0.01823434

Rmax = 100
alpha <- alfa <- 0.05

filename = paste0("save/SS_features_Mfix_Zipfs_05",".Rdat")
load(filename)
Mmax_all = ExpRes_list[[1]][,1]

mean_grid = seq(1e-8,1-1e-8,length.out = 1000)
UB_mat = matrix(-1, nrow = length(mean_grid), ncol = 1)

for(i in seq_along(mean_grid)){
  vmax = mean_grid[i] * (1 - mean_grid[i])
  v = vmax * (1 - 2*tol)
  ab = compute_ab_beta(mean_grid[i],vmax/2)
  UB_mat[i,1] = min(1, exp(compute_log_UBMarkov_FB(Rmax,ab[1],ab[2],n,Kn,M,alpha) ) )
}


## axis labels
ymax = (11/10) * max(UB_mat,Mmax_all); ymin = (10/11) * min(UB_mat,Mmax_all)
ylim_plot = c(ymin,ymax)
ypos = seq(ymin,ymax,length.out = 5)
ylabs = as.character(round(ypos*1e3,0))
xmax = max(mean_grid); xmin = min(mean_grid)
xpos = seq(xmin,xmax,length.out = 10)
xlabs = as.character(floor(xpos*100))
xlim_plot = c(xmin,xmax)

par( mfrow = c(1,1), mar = c(3.5,4.25,1,1), mgp=c(2.75,1,0), bty = "l", las = 1, cex.lab = 2 )
plot(0,0,  yaxt = "n", xaxt = "n",
     xlab = "", ylab = "1000 * bound",
     xlim = xlim_plot , ylim = ylim_plot, 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
mtext("mean * 100", side = 1, line = 2.5, cex = cex.axis)
points( x = mean_grid, y = UB_mat[,1], 
        type = "l", lwd = 5, col = "black" )
abline(h = Mmax_all, lty = 2, lwd = 0.5, col = "green")











