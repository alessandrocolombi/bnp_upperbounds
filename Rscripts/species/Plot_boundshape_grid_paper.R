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
save_img = TRUE
width = 22  # width of the output image
height = 10  # height of the output image
cex_plot = 1.8   # overall expansion factor for plot text/elements
cex_axis = 1.25  # size of axis tick labels
cex_lab  = 1.5   # size of y-axis titles
cex_xlab = 1.85  # size of x-axis titles only
lwd_curve = 4  # thickness of plotted curves
mfrow_plot = c(2, 3)  # plot layout: number of rows and columns
mar_plot = c(3.5, 4, 1, 0.6)  # margins: bottom, left, top, right
mgp_plot = c(1.8, 0.7, 0)  # axis layout: title line, tick-label line, axis line


# Options -----------------------------------------------------------------
n = 500
Rmax = 100
alpha <- alfa <- 0.05
M_max = 200
Kn_grid_pd = c(20, 100, 200, 400)
top_grid = data.frame(
  M = c(150, 200, 400, 500),
  Kn = c(50, 50, 200, 200)
)
curve_cols_top = hcl.colors(nrow(top_grid), palette = "viridis")
curve_cols_pd = hcl.colors(length(Kn_grid_pd), palette = "viridis")
top_legend = apply(top_grid, 1, function(x) {
  paste0("M=", x[["M"]], ", Kn=", x[["Kn"]])
})
pd_legend = paste0("Kn=", Kn_grid_pd)

fdp_grid_list = list(
  "3" = subset(
    expand.grid(M = c(500),
                Kn = c(50, 100, 150, 200)),
    Kn < M
  ),
  "6" = subset(
    expand.grid(M = c(500),
                Kn = c(50, 100, 150, 200)),
    Kn < M
  )
)

plot_dm_panel = function(ymax = 18) {
  gamma_grid = c(seq(0.001, 1, length.out = 100),
                 seq(1, 150, length.out = 200))
  y_ticks = seq(2, ymax, length.out = 5)
  plot(NA_real_, NA_real_, type = "n",
       xlim = c(0, 150), ylim = c(2, ymax),
       xlab = "", ylab = expression(q[M] == delta[M[0]]), yaxt = "n")
  axis(side = 2, at = y_ticks, labels = round(y_ticks, 0), cex.axis = cex_axis)
  for (ii in seq_len(nrow(top_grid))) {
    M = top_grid$M[ii]
    Kn = top_grid$Kn[ii]
    U_DM_grid = rep(NA_real_, length(gamma_grid))
    for (hh in seq_along(gamma_grid)) {
      gamma = gamma_grid[hh]
      ub_dm = exp(compute_log_UB_DirMulti(Rmax, gamma, M, Kn, n, alpha))
      U_DM_grid[hh] = min(ub_dm, 1)
    }
    lines(gamma_grid, 1000 * U_DM_grid, lwd = lwd_curve, col = curve_cols_top[ii])
  }
  legend("topright", legend = top_legend, col = curve_cols_top, lwd = lwd_curve,
         bty = "n", cex = 0.75)
  mtext(expression(gamma), side = 1, line = 1.5, cex = cex_xlab)
}

plot_fdp_panel = function(idx_lambda, ymax = 18) {
  Lambda_grid = c(50, 150, 200, 250, 300, 500)
  gamma_grid = c(seq(0.001, 1, length.out = 100),
                 seq(1, 150, length.out = 200))
  y_ticks = seq(2, ymax, length.out = 5)
  Lambda = Lambda_grid[idx_lambda]
  fdp_grid = fdp_grid_list[[as.character(idx_lambda)]]
  if (is.null(fdp_grid)) {
    stop("No FDP grid configured for idx_lambda = ", idx_lambda)
  }
  fdp_cols = hcl.colors(nrow(fdp_grid), palette = "viridis")
  fdp_legend = apply(fdp_grid, 1, function(x) {
    paste0("M=", x[["M"]], ", Kn=", x[["Kn"]])
  })
  plot(NA_real_, NA_real_, type = "n",
       ylim = c(2, ymax),
       xlim = range(gamma_grid),
       xlab = "",
       ylab = bquote(q[M] == Pois[1] * "(" * .(Lambda) * ")"),
       yaxt = "n")
  axis(side = 2, at = y_ticks, labels = round(y_ticks, 0), cex.axis = cex_axis)
  for (ii in seq_len(nrow(fdp_grid))) {
    M_val = fdp_grid$M[ii]
    Kn = fdp_grid$Kn[ii]
    ub_vals = rep(NA_real_, length(gamma_grid))
    for (hh in seq_along(gamma_grid)) {
      gamma = gamma_grid[hh]
      ub_fdp = exp(compute_log_UBMarkov_FD(Rmax, gamma, Lambda, Kn, n, alpha, M_val))
      ub_vals[hh] = min(ub_fdp, 1)
    }
    lines(gamma_grid, 1000 * ub_vals, lwd = lwd_curve, col = fdp_cols[ii])
  }
  legend("topright", legend = fdp_legend, col = fdp_cols, lwd = lwd_curve,
          bty = "n", cex = 0.7)
  mtext(expression(gamma), side = 1, line = 1.5, cex = cex_xlab)
}

plot_pd_panel = function(idx_sigma, ymax = 12.5) {
  sigma_grid = c(0, 0.25, 0.5, 0.7, 0.9, 0.99)
  theta_grid = c(seq(1e-10, 1, length.out = 1000),
                 seq(1, 1000, length.out = 1000))
  sigma = sigma_grid[idx_sigma]
  x_ticks = c(theta_grid[1], 350, 650, theta_grid[length(theta_grid)])
  x_tick_labels = c("0", "350", "650", "1000")
  y_ticks = seq(2, ymax, length.out = 5)

  plot(NA_real_, NA_real_, type = "n",
       ylim = c(2, ymax),
       xlim = range(theta_grid),
       xlab = "", xaxt = "n",
       ylab = bquote(sigma == .(sigma)), yaxt = "n")
  axis(side = 2, at = y_ticks, labels = round(y_ticks, 0), cex.axis = cex_axis)
  for (ii in seq_along(Kn_grid_pd)) {
    Kn = Kn_grid_pd[ii]
    ub_vals = rep(NA_real_, length(theta_grid))
    for (hh in seq_along(theta_grid)) {
      theta = theta_grid[hh]
      ub_pd = exp(compute_log_UBMarkov(Rmax, sigma, theta, Kn, n, alpha))
      ub_vals[hh] = min(ub_pd, 1)
    }
    lines(theta_grid, 1000 * ub_vals, lwd = lwd_curve, col = curve_cols_pd[ii])
  }
  if (idx_sigma == 1) {
    legend("topright", legend = pd_legend, col = curve_cols_pd, lwd = lwd_curve,
           bty = "n", cex = 0.75)
  }
  axis(side = 1, at = x_ticks, labels = x_tick_labels, cex.axis = cex_axis)
  mtext(expression(theta), side = 1, line = 1.5, cex = cex_xlab)
}

if (save_img) {
  pdf("img/UB_shape_grid_paper.pdf", width = width, height = height)
}

par(
  mfrow = mfrow_plot,
  mar = mar_plot,
  mgp = mgp_plot,
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

# # Row 2: same three PD panels as in UB_shape_PD.pdf
plot_pd_panel(1)
plot_pd_panel(3)
plot_pd_panel(5)

if (save_img) {
  dev.off()
}
