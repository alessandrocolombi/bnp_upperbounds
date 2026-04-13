# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc, wd_unicatt, wd_g100, wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd, "bnp_upperbounds/Rscripts/species")
setwd(wd)

# Functions ---------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))

source("../../R/Rfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")

# Plot options ------------------------------------------------------------
save_img = FALSE
width = 22
height = 10
cex_plot = 1.8  # overall expansion factor for plot text/elements
cex_axis = 1.25  # size of axis tick labels
cex_lab  = 1.5  # size of axis titles
cex_xlab = 1.85  # size of x-axis titles only
lwd_curve = 4
curve_cols = c("red2", "darkorange", "deepskyblue3", "darkviolet")

# Options -----------------------------------------------------------------
n = 500
Rmax = 100
alpha <- alfa <- 0.05
M_max = 200
M_grid_top = c(200, 400)
Kn_grid_top = c(50, 100)
Kn_grid_pd = c(20, 100, 200, 400)
top_grid = expand.grid(M = M_grid_top, Kn = Kn_grid_top)
top_legend = apply(top_grid, 1, function(x) {
  paste0("M=", x[["M"]], ", Kn=", x[["Kn"]])
})
pd_legend = paste0("Kn=", Kn_grid_pd)

plot_dm_panel = function() {
  gamma_grid = c(seq(0.001, 1, length.out = 100),
                 seq(1, 150, length.out = 200))
  plot(NA_real_, NA_real_, type = "n",
       xlim = c(0, 150), ylim = c(2, 14),
       xlab = "", ylab = expression(q[M] == delta[M[0]]))
  for (ii in seq_len(nrow(top_grid))) {
    M = top_grid$M[ii]
    Kn = top_grid$Kn[ii]
    U_DM_grid = rep(NA_real_, length(gamma_grid))
    for (hh in seq_along(gamma_grid)) {
      gamma = gamma_grid[hh]
      ub_dm = exp(compute_log_UB_DirMulti(Rmax, gamma, M, Kn, n, alpha))
      U_DM_grid[hh] = min(ub_dm, 1)
    }
    lines(gamma_grid, 1000 * U_DM_grid, lwd = lwd_curve, col = curve_cols[ii])
  }
  legend("topright", legend = top_legend, col = curve_cols, lwd = lwd_curve,
         bty = "n", cex = 0.75)
  mtext(expression(gamma), side = 1, line = 1.5, cex = cex_xlab)
}

plot_fdp_panel = function(idx_lambda) {
  Lambda_grid = c(50, 150, 200, 250, 300, 500)
  gamma_grid = c(seq(0.001, 1, length.out = 100),
                 seq(1, 150, length.out = 200))
  Lambda = Lambda_grid[idx_lambda]
  plot(NA_real_, NA_real_, type = "n",
       ylim = c(2, 40),
       xlim = range(gamma_grid),
       xlab = "",
       ylab = bquote(q[M] == Pois[1] * "(" * .(Lambda) * ")"))
  for (ii in seq_len(nrow(top_grid))) {
    M_val = top_grid$M[ii]
    Kn = top_grid$Kn[ii]
    ub_vals = rep(NA_real_, length(gamma_grid))
    for (hh in seq_along(gamma_grid)) {
      gamma = gamma_grid[hh]
      ub_fdp = exp(compute_log_UBMarkov_FD(Rmax, gamma, Lambda, Kn, n, alpha, M_val))
      ub_vals[hh] = min(ub_fdp, 1)
    }
    lines(gamma_grid, 1000 * ub_vals, lwd = lwd_curve, col = curve_cols[ii])
  }
  mtext(expression(gamma), side = 1, line = 1.5, cex = cex_xlab)
}

plot_pd_panel = function(idx_sigma) {
  sigma_grid = c(0, 0.25, 0.5, 0.7, 0.9, 0.99)
  theta_grid = c(seq(1e-10, 1, length.out = 1000),
                 seq(1, 1000, length.out = 1000))
  sigma = sigma_grid[idx_sigma]
  x_ticks = c(theta_grid[1], 350, 650, theta_grid[length(theta_grid)])
  x_tick_labels = c("0", "350", "650", "1000")

  plot(NA_real_, NA_real_, type = "n",
       ylim = c(2, 12.5),
       xlim = range(theta_grid),
       xlab = "", xaxt = "n",
       ylab = bquote(sigma == .(sigma)))
  for (ii in seq_along(Kn_grid_pd)) {
    Kn = Kn_grid_pd[ii]
    ub_vals = rep(NA_real_, length(theta_grid))
    for (hh in seq_along(theta_grid)) {
      theta = theta_grid[hh]
      ub_pd = exp(compute_log_UBMarkov(Rmax, sigma, theta, Kn, n, alpha))
      ub_vals[hh] = min(ub_pd, 1)
    }
    lines(theta_grid, 1000 * ub_vals, lwd = lwd_curve, col = curve_cols[ii])
  }
  if (idx_sigma == 1) {
    legend("topright", legend = pd_legend, col = curve_cols, lwd = lwd_curve,
           bty = "n", cex = 0.75)
  }
  axis(side = 1, at = x_ticks, labels = x_tick_labels, cex.axis = cex_axis)
  mtext(expression(theta), side = 1, line = 1.5, cex = cex_xlab)
}

if (save_img) {
  pdf("img/UB_shape_grid_paper.pdf", width = width, height = height)
}

par(
  mfrow = c(2, 3),
  mar = c(3.5, 6.2, 1, 1),
  mgp = c(1.8, 0.7, 0),
  bty = "l",
  las = 1,
  cex = cex_plot,
  cex.axis = cex_axis,
  cex.lab = cex_lab
)

# Row 1: DM, FDP position 2, FDP position 3
plot_dm_panel()
plot_fdp_panel(3)
plot_fdp_panel(6)

# Row 2: same three PD panels as in UB_shape_PD.pdf
plot_pd_panel(1)
plot_pd_panel(3)
plot_pd_panel(5)

if (save_img) {
  dev.off()
}
