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

if (save_img) {
  pdf("img/UB_shape_features_paper.pdf", width = width, height = height)
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

# Finite Beta plot --------------------------------------------------------

Rmax = 100
n = 5000
M = 5000
Mstar_grid = c(1, 50, 500)
alpha <- alfa <- 0.05
tol = sqrt(.Machine$double.eps) 

mean_grid = seq(1e-5, 1 - 1e-5, length.out = 1000)
kappa_grid = c(2, 10, 50, 100)
mycol_ub = hcl.colors(n = length(kappa_grid), palette = "viridis")
xmax = max(mean_grid)
xmin = min(mean_grid)
xpos = seq(xmin, xmax, length.out = 6)
xlabs = as.character(round(xpos, 2))
xlim_plot = c(xmin, xmax)
UB_list = vector("list", length(Mstar_grid))

for (ii in seq_along(Mstar_grid)) {
  Mstar = Mstar_grid[ii]
  Kn = M - Mstar
  UB_mat = matrix(-1, nrow = length(mean_grid), ncol = length(kappa_grid))

  for (i in seq_along(mean_grid)) {
    vmax = mean_grid[i] * (1 - mean_grid[i])
    for (j in seq_along(kappa_grid)) {
      v = vmax / kappa_grid[j]
      ab = compute_ab_beta(mean_grid[i], v)
      UB_mat[i, j] = min(1, exp(compute_log_UBMarkov_FB(Rmax, ab[1], ab[2], n, Kn, M, alpha)))
    }
  }
  UB_list[[ii]] = UB_mat
}

global_ymax = (11 / 10) * max(unlist(UB_list))
global_ymin = (10 / 11) * min(unlist(UB_list))
ylim_plot = c(global_ymin, global_ymax)
ypos = seq(global_ymin, global_ymax, length.out = 5)
ylabs = as.character(round(1000 * ypos, 0))



for (ii in seq_along(Mstar_grid)) {
  Mstar = Mstar_grid[ii]
  UB_mat = UB_list[[ii]]

  plot(0, 0, yaxt = "n", xaxt = "n",
       xlab = "", ylab = bquote(M^{"*"} == .(Mstar)),
       xlim = xlim_plot, ylim = ylim_plot,
       main = " ",
       type = "n")
  grid(lty = 1, lwd = 1, col = "gray90")
  axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex_axis)
  axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex_axis)
  mtext(expression(omega), side = 1, line = 2.2, cex = cex_xlab)
  for (j in seq_along(kappa_grid)) {
    points(x = mean_grid, y = UB_mat[, j],
           type = "l", lwd = lwd_curve, col = mycol_ub[j])
  }
  legend("topleft",
         legend = sapply(kappa_grid, function(x) as.expression(bquote(eta == .(x)))),
         col = mycol_ub, lwd = lwd_curve,
         bty = "n", cex = 0.8)
}



# IBP plot --------------------------------------------------------

Rmax = 100
n = 100
alfa = 0.05
gamma_grid = c(1, 100, 5000)
alpha_grid = c(0, 0.25, 0.5, 0.9)
c_grid = c(seq(1e-5, 1, length.out = 500),
           seq(1, 100, length.out = 500))
mycol_ibp = hcl.colors(n = length(alpha_grid), palette = "viridis")
xmax = max(c_grid)
xmin = min(c_grid)
xpos = c(0, 1, 25, 50, 75, 100)
xlabs = as.character(xpos)
xlim_plot = c(xmin, xmax)
UB_list = vector("list", length(gamma_grid))

for (ii in seq_along(gamma_grid)) {
  gamma = gamma_grid[ii]
  UB_mat = matrix(-1, nrow = length(c_grid), ncol = length(alpha_grid))

  for (i in seq_along(c_grid)) {
    c_val = c_grid[i]
    for (j in seq_along(alpha_grid)) {
      alpha_val = alpha_grid[j]
      UB_mat[i, j] = min(1, exp(compute_log_UBMarkov_BeBePois(Rmax, alpha_val, c_val, gamma, n, alfa)))
    }
  }
  UB_list[[ii]] = UB_mat
}

global_ymax = (11 / 10) * max(unlist(UB_list))
global_ymin = (10 / 11) * min(unlist(UB_list))
ylim_plot = c(global_ymin, global_ymax)
ypos = seq(global_ymin, global_ymax, length.out = 5)
ylabs = as.character(round(1000 * ypos, 0))


for (ii in seq_along(gamma_grid)) {
  gamma = gamma_grid[ii]
  UB_mat = UB_list[[ii]]

  plot(0, 0, yaxt = "n", xaxt = "n",
       xlab = "", ylab = bquote(gamma == .(gamma)),
       xlim = xlim_plot, ylim = ylim_plot,
       main = " ",
       type = "n")
  grid(lty = 1, lwd = 1, col = "gray90")
  axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex_axis)
  axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex_axis)
  mtext("c", side = 1, line = 2.2, cex = cex_xlab)
  for (j in seq_along(alpha_grid)) {
    points(x = c_grid, y = UB_mat[, j],
           type = "l", lwd = lwd_curve, col = mycol_ibp[j])
  }
  legend("topright",
         legend = sapply(alpha_grid, function(x) as.expression(bquote(sigma == .(x)))),
         col = mycol_ibp, lwd = lwd_curve,
         bty = "n", cex = 0.8)
}


if (save_img) {
  dev.off()
}


