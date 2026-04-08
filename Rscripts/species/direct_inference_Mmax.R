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
source("../../R/PFFAfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")


# Custom functions --------------------------------------------------------

run_single_rep = function(model, n, Kn, params, B, seed = 123231,
                          alpha = 0.05, chunk_size = 5000, Natoms_max=1e5)
{
  library(matrixStats)
  Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
  set.seed(seed)
  
  if(model == "FDP"){
    gamma = params[1]
    Mstar = as.integer(params[2])
    n_chunks = ceiling(B / chunk_size)
    x_all = numeric(B)
    beta_MC = rbeta(n = B,shape1 = Mstar * gamma,shape2 = n + Kn * gamma)
    start = 1
    seeds = sample.int(999999, size = n_chunks)
    for(cc in seq_len(n_chunks)){
      end = min(cc * chunk_size, B)
      m = end - start + 1
      
      # G = r_gamma_mat(m, Mstar, gamma, 1.0, seeds[cc])
      G = matrix(rgamma(n = m*Mstar,shape=gamma,rate=1),nrow = m, ncol = Mstar)
      vmax = rowMaxs(G) / rowSums(G)
      
      x_all[start:end] = beta_MC[start:end] * vmax
      start = end + 1
    }
    return(unname(quantile(x_all, 1 - alpha, na.rm = TRUE)))
  }
  
  if(model == "PYP"){
    sigma = params[1]; theta = params[2]
    x_all = numeric(B)
    beta_MC = rbeta(n = B, shape1 = theta + Kn*sigma,
                    shape2 = n - Kn*sigma)
    Vmax = r_SB_max(Nrep=B,Natoms_max=Natoms_max,
                    alpha=sigma,theta=theta,
                    seed=seed)$pmax
    x_all = Vmax*beta_MC
    return(unname(quantile(x_all, 1 - alpha, na.rm = TRUE)))
  }
  
  stop("Model must either be FDP or PYP")
}
# Plot options ------------------------------------------------------------------
save_img = TRUE
width = 12; height = 6
cex.labels <- cex.lab <- 2
cex.axis <- 2
cex.legend <- 1.5
mycol = c("darkgreen","darkorange","darkred","darkblue","lightblue")
mycol2 = c("black","lightblue")
lgd_names = c("Unbounded","Bounded","IBP","MBP","FB")

# FDP case ----------------------------------------------------------------
alpha = 0.05
## Approximate reasoning -----------------------------------------------
n = 1e4
thresholds = c(0.218,0.069,0.0218)
omega_grid = seq(1e-3,5,length.out = 1000)
qfq = rep(-1,length(omega_grid))
for(i in seq_along(omega_grid)){
  omega = omega_grid[i]
  q = qbeta(p=1-alpha, shape1 = omega, shape2 = n)
  fq = dbeta(x = q, shape1 = omega, shape2 = n )
  qfq[i] = q*fq
}

par( mfrow = c(1,1), mar = c(4,7,1,2), mgp=c(5,1,0), bty = "l", 
     las = 1, cex.lab = cex.lab, cex.axis = cex.axis )
plot(x = omega_grid, y = qfq, type = "l", lwd = 3,
     xlab = "", ylab = "target")
abline(h = thresholds, col = "red", lty = 2, lwd = 2)
mtext(expression(omega), side = 1, line = 2.5, cex = cex.axis)


omega_target = sapply(thresholds, function(t) {
  idx=which(qfq<t);omega_grid[max(idx)]
})
omega_target

## Generate data -----------------------------------------------------------
seed = 132131
set.seed(seed)
n = 1000
M = 1000 #<---
gamma = omega_target[1]/M #<---
gamma = 0.005
Kn = 50
Mstar = M-Kn

cat("\n",gamma," || ",Kn," || ", Mstar )


## Compute true quantile -----------------------------------------------
library(matrixStats)
tictoc::tic()
B4true = 1e6
chunk_size = 5000
n_chunks = ceiling(B4true / chunk_size)

x_all = numeric(B4true)
beta_MC = rbeta(n = B4true, shape1 = Mstar * gamma, shape2 = n + Kn * gamma)

start = 1
seeds = sample(1:999999,size = n_chunks)
for(cc in seq_len(n_chunks)){
  end = min(cc * chunk_size, B4true)
  m = end - start + 1
  G = r_gamma_mat(m,Mstar,gamma,1.0,seeds[cc])
  vmax = rowMaxs(G) / rowSums(G)
  x_all[start:end] = beta_MC[start:end] * vmax
  start = end + 1
}
tictoc::toc()

q_true = unname(quantile(x_all, 1 - alpha, names = FALSE))
dens_obj = density(x_all, n = 2^14)
plot(x = dens_obj$x, y = dens_obj$y); abline(v = q_true, lty = 2, col = "red")
fhat_kde = approx(x = dens_obj$x, y = dens_obj$y, xout = q_true, rule = 2)$y

cat("\n",q_true," || ",fhat_kde," || ",fhat_kde*q_true, "\n" )


## Simulation options --------------------------------------------------
model = "FDP"
num_cores = 3 #<---
Nrep = 100 #<---
Bgrid = c(1e1,1e2,1e3,1e4) #<---
LBgrid = length(Bgrid)
qest_all = matrix(-1,nrow = Nrep, ncol = LBgrid)

ii = 2
for(ii in seq_along(Bgrid)){
  B = Bgrid[ii]
  seeds = sample(1:999999, size = Nrep)
  cat("\n","B = ",B,"\n")
  cluster <- makeCluster(num_cores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)
  clusterExport(cluster, list(), envir = environment())
  res_list = parLapply( cl = cluster, x = seeds,
                        fun = run_single_rep,
                        model=model, n=n, Kn=Kn, 
                        params=c(gamma,Mstar), B=B, alpha=alpha )
  stopCluster(cluster)
  
  # Compute relevant summaries
  qest_all[,ii] = unlist(res_list)

}

summaries = data.frame(
  B      = Bgrid,
  mean   = colMeans(qest_all),
  median = apply(qest_all, 2, median),
  sd     = apply(qest_all, 2, sd),
  var    = apply(qest_all, 2, var),
  bias   = colMeans(qest_all) - q_true,
  rmse   = sqrt(colMeans((qest_all - q_true)^2)),
  rel_rmse   = sqrt(colMeans((qest_all - q_true)^2))/ q_true,
  p05    = apply(qest_all, 2, quantile, probs = 0.05, names = FALSE),
  p95    = apply(qest_all, 2, quantile, probs = 0.95, names = FALSE)
)

print(summaries)


## Plot Mario ------------------------------------------------------------------
xlim_plot = c( min(Bgrid), max(Bgrid) )
ylim_plot = range(c(summaries$p05, summaries$p95, q_true))
ypos = seq(ylim_plot[1],ylim_plot[2], length.out = 5)
ylabs = round(ypos*1e3,1)
xpos  = Bgrid
xlabs = format(xpos, scientific = FALSE, trim = TRUE)


if(save_img)
  pdf("img/DirectMmax_FDP_EstQ.pdf", width=width, height=height)
par( mfrow = c(1,1), mar = c(4,5,1,2), mgp=c(3.5,1,0), bty = "l", 
     las = 1, cex.lab = cex.lab, cex.axis = cex.axis )
plot( x = Bgrid, y = summaries$median, log = "x",
      type = "b", pch = 16, lwd = 3,
      xaxt = "n", yaxt = "n",
      xlab = "", ylab = "Est. quantile",
      xlim = xlim_plot , ylim = ylim_plot )
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
mtext("B", side = 1, line = 2.5, cex = cex.axis)
polygon( c(Bgrid, rev(Bgrid)),
         c(summaries$p05, rev(summaries$p95)),
         col = scales::alpha("lightblue", 0.25),
         border = NA) # plot in-sample bands
if(save_img)
  dev.off()


## Plot RMSE ---------------------------------------------------------------
xlim_plot = c( min(Bgrid), max(Bgrid) )
ylim_plot = c(0.05,max(summaries$rel_rmse))#range(c(summaries$rel_rmse))
ypos = summaries$rel_rmse
ylabs = format(round(ypos*1e2,1), scientific = FALSE, trim = FALSE)
xpos  = Bgrid
xlabs = format(xpos, scientific = FALSE, trim = TRUE)

if(save_img)
  pdf("img/DirectMmax_FDP_RelRMSE.pdf", width=width, height=height)
par( mfrow = c(1,1), mar = c(4,6,1,2), mgp=c(4.25,1,0), bty = "l", 
     las = 1, cex.lab = cex.lab, cex.axis = cex.axis )
plot( x = Bgrid, y = summaries$rel_rmse, log = "xy",
      type = "b", lwd = 3, pch = 16,
      xaxt = "n", yaxt = "n",
      xlab = "", ylab = "Rel. RMSE" )
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
mtext("B", side = 1, line = 2.5, cex = cex.axis)
abline(h = 0.01, lty = 2)
if(save_img)
  dev.off()

# sqrt(B) reference line
c0 = summaries$rel_rmse[1] * sqrt(Bgrid[1])
ref_line = c0 / sqrt(Bgrid)  
lines(Bgrid, ref_line, lty = 3, lwd = 2)



# PYP case ----------------------------------------------------------------
alpha = 0.05
seed = 132131
set.seed(seed)
## Generate data -----------------------------------------------------------
seed = 132131
set.seed(seed)
n = 1000
sigma = 0;Kn = 50
theta = 3

SB = r_SB(Nrep = 10, Natoms = 10000, alpha = sigma, theta = theta, seed = 3121)
round(SB[,1:5],4)

## Compute true quantile -----------------------------------------------
Natoms_max = 10000
Vmax_all = r_SB_max(Nrep=10,Natoms_max=Natoms_max,alpha=sigma,theta=theta,seed=seed)
Vmax_all$converged
Vmax_all$n_atoms_used

tictoc::tic()
B4true = 1e6
beta_MC = rbeta(n = B4true, shape1 = theta + Kn*sigma,
                shape2 = n - Kn*sigma)
Vmax = r_SB_max(Nrep=B4true,Natoms_max=Natoms_max,
                alpha=sigma,theta=theta,
                seed=seed)$pmax
x_all = Vmax*beta_MC
tictoc::toc()

q_true = unname(quantile(x_all, 1 - alpha, names = FALSE))
dens_obj = density(x_all, n = 2^14)
plot(x = dens_obj$x, y = dens_obj$y); abline(v = q_true, lty = 2, col = "red")
fhat_kde = approx(x = dens_obj$x, y = dens_obj$y, xout = q_true, rule = 2)$y

cat("\n",q_true," || ",fhat_kde," || ",fhat_kde*q_true, "\n" )


## Simulation options --------------------------------------------------
model = "PYP"
num_cores = 3 #<---
Nrep = 10 #<---
Bgrid = c(1e1,1e2,1e3,1e4) #<---
LBgrid = length(Bgrid)
qest_pyp_all = matrix(-1,nrow = Nrep, ncol = LBgrid)

ii = 2
for(ii in seq_along(Bgrid)){
  B = Bgrid[ii]
  seeds = sample(1:999999, size = Nrep)
  cat("\n","B = ",B,"\n")
  cluster <- makeCluster(num_cores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)
  clusterExport(cluster, list("r_SB_max"), envir = environment())
  res_list = parLapply( cl = cluster, x = seeds,
                        fun = run_single_rep,
                        model=model, n=n, Kn=Kn, 
                        params=c(sigma,theta), B=B, alpha=alpha, Natoms_max=Natoms_max )
  stopCluster(cluster)
  
  # Compute relevant summaries
  qest_pyp_all[,ii] = unlist(res_list)
  
}

summaries = data.frame(
  B      = Bgrid,
  mean   = colMeans(qest_pyp_all),
  median = apply(qest_pyp_all, 2, median),
  sd     = apply(qest_pyp_all, 2, sd),
  var    = apply(qest_pyp_all, 2, var),
  bias   = colMeans(qest_pyp_all) - q_true,
  rmse   = sqrt(colMeans((qest_pyp_all - q_true)^2)),
  rel_rmse   = sqrt(colMeans((qest_pyp_all - q_true)^2))/ q_true,
  p05    = apply(qest_pyp_all, 2, quantile, probs = 0.05, names = FALSE),
  p95    = apply(qest_pyp_all, 2, quantile, probs = 0.95, names = FALSE)
)

print(summaries)


## Plot Mario ------------------------------------------------------------------
xlim_plot = c( min(Bgrid), max(Bgrid) )
ylim_plot = range(c(summaries$p05, summaries$p95, q_true))
ypos = seq(ylim_plot[1],ylim_plot[2], length.out = 5)
ylabs = round(ypos*1e3,1)
xpos  = Bgrid
xlabs = format(xpos, scientific = FALSE, trim = TRUE)


if(save_img)
  pdf("img/DirectMmax_PYP_EstQ.pdf", width=width, height=height)
par( mfrow = c(1,1), mar = c(4,5,1,2), mgp=c(3.5,1,0), bty = "l", 
     las = 1, cex.lab = cex.lab, cex.axis = cex.axis )
plot( x = Bgrid, y = summaries$median, log = "x",
      type = "b", pch = 16, lwd = 3,
      xaxt = "n", yaxt = "n",
      xlab = "", ylab = "Est. quantile",
      xlim = xlim_plot , ylim = ylim_plot )
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
mtext("B", side = 1, line = 2.5, cex = cex.axis)
polygon( c(Bgrid, rev(Bgrid)),
         c(summaries$p05, rev(summaries$p95)),
         col = scales::alpha("lightblue", 0.25),
         border = NA) # plot in-sample bands
if(save_img)
  dev.off()


## Plot RMSE ---------------------------------------------------------------
xlim_plot = c( min(Bgrid), max(Bgrid) )
ylim_plot = c(0.05,max(summaries$rel_rmse))#range(c(summaries$rel_rmse))
ypos = summaries$rel_rmse
ylabs = format(round(ypos*1e2,1), scientific = FALSE, trim = FALSE)
xpos  = Bgrid
xlabs = format(xpos, scientific = FALSE, trim = TRUE)

if(save_img)
  pdf("img/DirectMmax_PYP_RelRMSE.pdf", width=width, height=height)
par( mfrow = c(1,1), mar = c(4,6,1,2), mgp=c(4.25,1,0), bty = "l", 
     las = 1, cex.lab = cex.lab, cex.axis = cex.axis )
plot( x = Bgrid, y = summaries$rel_rmse, log = "xy",
      type = "b", lwd = 3, pch = 16,
      xaxt = "n", yaxt = "n",
      xlab = "", ylab = "Rel. RMSE" )
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
mtext("B", side = 1, line = 2.5, cex = cex.axis)
abline(h = 0.01, lty = 2)
if(save_img)
  dev.off()

# sqrt(B) reference line
# c0 = summaries$rel_rmse[1] * sqrt(Bgrid[1])
# ref_line = c0 / sqrt(Bgrid)  
# lines(Bgrid, ref_line, lty = 3, lwd = 2)













