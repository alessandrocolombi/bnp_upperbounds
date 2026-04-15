# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc, wd_unicatt, wd_g100, wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd, "bnp_upperbounds/Rscripts/features")
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


## 1) Sanity check on synthetic Bernoulli-Beta data -----------------------
set.seed(42)
n = 2000
M = 100
a_true = 0.5
b_true = 70
p_true = rbeta(M, shape1 = a_true, shape2 = b_true)
n_i = rbinom(n = M, size = n, prob = p_true)
Kn = sum(n_i > 0)
hy_FB = make_hy_FB(n = n, M = M, Kn = Kn, mu_kappa = mu_kappa, var_kappa = 10 * max(M - Kn, 1))
a_kappa = hy_FB[3]
b_kappa = hy_FB[4]

kappa_grid = seq(1e-5, 1000, length.out = 10000)
plot(x = kappa_grid, y = dgamma(x = kappa_grid, shape = a_kappa, rate = b_kappa),
     xlab = "kappa", ylab = "")

fit_synth = fit_FB_map(n = n, M = M, n_i = n_i, hy = hy_FB)

cat("\n=== Synthetic Bernoulli-Beta sanity check ===\n")
cat("True a =", a_true, "| True b =", b_true, "\n")
cat("True mu =", a_true/(a_true+b_true), "| True b =", a_true+b_true+1, "\n")
cat("MAP a  =", round(fit_synth$a_map, 4), "| MAP b  =", round(fit_synth$b_map, 4), "\n")
cat("MAP mu =", round(fit_synth$mu_map, 4), "| MAP kappa =", round(fit_synth$kappa_map, 4), "\n")
cat("Observed Kn =", fit_synth$Kn, "\n")


## 2) Zipfs data with MAP fit for FB only ---------------------------------
set.seed(123)
n = 2000
Rmax = 100
alpha <- alfa <- 0.05
name = "Zipfs"
params = 0.85
M_grid = c(1100, 1600, 2100, 3600)
Nrep = 50


mu_kappa = 200
Zipfs_FB_MAP = lapply(M_grid, function(M) {
  ptrue = sim_features_generic(name, M, params)
  ptrue = sort(ptrue, decreasing = TRUE)

  rep_res = lapply(seq_len(Nrep), function(rep_idx) {
    n_i = rbinom(n = length(ptrue), size = n, prob = ptrue)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)

    if (Kn == M) {
      Mmax = 0
    } else {
      idx_unobs = which(n_i == 0)
      Mmax = max(ptrue[idx_unobs])
    }

    hy_FB = make_hy_FB( n = n, M = M, Kn = Kn, mu_kappa = mu_kappa, var_kappa=10*max(1,M-Kn) )
    fit_map = fit_FB_map(n = n, M = M, n_i = n_i, hy = hy_FB)
    a_map = fit_map$a_map
    b_map = fit_map$b_map

    ubFB = exp(compute_log_UBMarkov_FB(Rmax, a_map, b_map, n, Kn, M, alpha))
    ubFB = min(ubFB, 1)

    data.frame(
      rep = rep_idx,
      M = M,
      Kn = Kn,
      Mstar = M-Kn,
      Mmax = Mmax,
      a_map = a_map,
      b_map = b_map,
      mu_map = fit_map$mu_map,
      kappa_map = fit_map$kappa_map,
      ubFB = ubFB,
      covered = as.integer(Mmax <= ubFB),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rep_res)
})
Zipfs_FB_MAP = do.call(rbind, Zipfs_FB_MAP)
Zipfs_FB_MAP[1:10,]
Zipfs_FB_MAP[11:21,]
Zipfs_FB_MAP[21:30,]
Zipfs_FB_MAP[31:40,]

 
Zipfs_FB_MAP_cov = do.call(rbind, lapply(split(Zipfs_FB_MAP, Zipfs_FB_MAP$M), function(x) {
  data.frame(
    M = x$M[1],
    coverage = mean(x$covered),
    mean_Kn = mean(x$Kn),
    mean_Mmax = mean(x$Mmax),
    mean_ubFB = mean(x$ubFB),
    mean_a_map = mean(x$a_map),
    mean_b_map = mean(x$b_map),
    stringsAsFactors = FALSE
  )
}))
rownames(Zipfs_FB_MAP_cov) = NULL

cat("\n=== Zipfs(0.85) FB-MAP fit: coverage over", Nrep, "replications ===\n")
print(Zipfs_FB_MAP_cov)


## 3) Zipfs likelihood slices in mu for fixed kappa ----------------------
set.seed(3231)
n = 2000
name = "Zipfs"
params = 0.85
M_like = 1100 # <--- modify here

ptrue = sim_features_generic(name, M_like, params)
ptrue = sort(ptrue, decreasing = TRUE)
n_i = rbinom(n = length(ptrue), size = n, prob = ptrue)
Kn = sum(n_i > 0); Kn

mu_grid = seq(1e-5, 0.2, length.out = 1000)
kappa_grid_plot = seq(10, 500, length.out = 6)
llik_mat = matrix(NA_real_, nrow = length(mu_grid), ncol = length(kappa_grid_plot))
mean_mu_grid = rep(NA_real_, length(kappa_grid_plot))

for (j in seq_along(kappa_grid_plot)) {
  kappa = kappa_grid_plot[j]
  mean_mu_grid[j] = find_mean_mu_FB(n = n, Kn = Kn, kappa = kappa, M = M_like)
  for (i in seq_along(mu_grid)) {
    mu = mu_grid[i]
    a = mu * (kappa - 1)
    b = (kappa - 1) * (1 - mu)
    llik_mat[i, j] = -llik_FB(x = c(a, b), n = n, Kn = Kn, data_obs = n_i, M = M_like)
  }
}

par(mfrow = c(2, 3), mar = c(3.5, 4.5, 2, 1), mgp = c(2.5, 0.8, 0), bty = "l", las = 1)
for (j in seq_along(kappa_grid_plot)) {
  yy = llik_mat[, j]
  ymax = max(yy, na.rm = TRUE)
  ymin = min(yy, na.rm = TRUE)
  ypos = seq(ymin, ymax, length.out = 5)
  ylabs = formatC(ypos, format = "f", digits = 1)
  xpos = seq(0, 1, length.out = 6)
  xlabs = formatC(xpos, format = "f", digits = 1)

  plot(0, 0, type = "n",
       xlim = c(0, max(mu_grid)), ylim = c(ymin, ymax),
       xlab = "", ylab = "log-likelihood",
       xaxt = "n", yaxt = "n",
       main = paste0("kappa = ", round(kappa_grid_plot[j], 1)))
  grid(lty = 1, lwd = 1, col = "gray90")
  axis(side = 1, at = xpos, labels = xlabs)
  axis(side = 2, at = ypos, labels = ylabs)
  mtext("mu", side = 1, line = 2.2)
  points(mu_grid, yy, type = "l", lwd = 3, col = "darkblue")
  abline(v = mean_mu_grid[j], lty = 2, lwd = 2, col = "red")
}







ExpKn_zipfs = function(M,n,s){
  w = sapply(2:(M+1),function(j) j^(-s))
  M - sum((1-w)^n)
}

ExpKn_zipfs(M=1100, n=2000, s=0.85)
ExpKn_zipfs(M=1600, n=2000, s=0.85)
ExpKn_zipfs(M=2100, n=2000, s=0.85)
ExpKn_zipfs(M=3600, n=2000, s=0.85)

ExpKn_zipfs(M=5000, n=2000, s=0.85)
ExpKn_zipfs(M=5000, n=7500, s=0.85)
ExpKn_zipfs(M=5000, n=9000, s=0.85)
ExpKn_zipfs(M=5000, n=10000, s=0.85)


ExpKn_zipfs(M=1100, n=2000, s=1.02)
ExpKn_zipfs(M=1600, n=2000, s=1.02)
ExpKn_zipfs(M=2100, n=2000, s=1.02)
ExpKn_zipfs(M=3600, n=2000, s=1.02)



c = c(2,1000,5000)
(1-(1/c)^n)
