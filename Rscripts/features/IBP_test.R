# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/CriminalData/")
setwd(wd)

# Functions ---------------------------------------------------------------
source("../../R/Rfunctions.R")
source("../../R/PFFAfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")

# Funzioni ----------------------------------------------------------------
base_rnd <- function(n) runif(n)
# Common params -----------------------------------------------------------
seed = 42
set.seed(seed)
alfa <- alpha <- 0.05
n = 100

# Sim. data from model -------------------------------------------------------
BB = 50

sigma_true = 0.25; c_true = 10; gamma_true = 10
xmax = 100

seeds = sample(1:999999, size = BB)
sim_3IBP_all = lapply(seeds, function(seed) fk_stable_beta_process( c = c_true, sigma = sigma_true, gamma = gamma_true,
                                                                    base_rnd = base_rnd,
                                                                    eps = 1e-6, seed = seed)$w )
sim_3IBP = lapply(sim_3IBP_all, function(x) x[1:xmax])
mat_3IBP = do.call(cbind, sim_3IBP)
qnt_3IBP = apply(mat_3IBP, 1, quantile, probs = c(0.025,0.5,0.975))

## Plot data
ymax = max(qnt_3IBP) * 1.02
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_3IBP[1,1:xmax], rev(qnt_3IBP[3,1:xmax])),
         col = "grey85",
         border = NA) # plot in-sample bands

# Gen. all datasets
lengths_vec <- lengths(sim_3IBP_all)   # L_i
Mmax <- max(lengths_vec)
row_idx <- rep(seq_along(lengths_vec), lengths_vec)   # row for each prob
col_idx <- sequence(lengths_vec)                       # column index within row
probs_all <- unlist(sim_3IBP_all)

# Draw all required binomials at once
draws <- rbinom(length(probs_all), size = n, prob = probs_all)

# Fill matrix (initialize with 0 or NA as you prefer)
data <- matrix(0L, nrow = length(lengths_vec), ncol = Mmax)
data[cbind(row_idx, col_idx)] <- draws
dim(data)

## Fit and check Ale -----------------------------------------------------------
# All datasets
start_params <- c(alpha = 0.1, gamma= 1, u = 1)
res_vec = apply(data, 1, function(counts){
  counts = counts[which(counts > 0)]
  Kn = length(counts)
  fit <- tryCatch(optim(par = start_params, fn = llik_PP3Parm, 
                        method = "L-BFGS-B",
                        n = n, Kn = Kn, data_obs = counts,
                        lower = c(1e-16, 1e-16, 1e-16), 
                        upper = c(1-1e-10, Inf, Inf)) ,
                  error=function(e) NA)
  list("sigma_mle" = fit$par[1],
       "c_mle" = fit$par[3] - fit$par[1],
       "gamma_mle" = fit$par[2],
       "fval" = fit$value)
})

res_mat = matrix(0,nrow = BB, ncol = 0)
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$sigma_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$c_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$gamma_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$fval) )
colnames(res_mat) <- names(res_vec[[1]])

colMeans(res_mat)

fitcheck_all = apply(res_mat, 1, function(x) {
  fk_stable_beta_process( c = x[2], sigma = x[1], 
                          gamma = x[3], base_rnd = base_rnd,
                          eps = 1e-6, seed = 12323123)$w 
})
fitcheck = lapply(fitcheck_all, function(x) x[1:xmax])
fitcheck = do.call(cbind, fitcheck); fitcheck[is.na(fitcheck)] = 0
fitcheck_qnt = apply(fitcheck, 1, quantile, probs = c(0.025,0.5,0.975))

## Plot data
ymax = max(max(qnt_3IBP),max(fitcheck_qnt)) * 1.02
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1)
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_3IBP[1,1:xmax], rev(qnt_3IBP[3,1:xmax])),
         col = "grey85",
         border = NA) # plot in-sample bands
polygon( c(1:xmax, rev(1:xmax)),
         c(fitcheck_qnt[1,1:xmax], rev(fitcheck_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands



## Fit and check -----------------------------------------------------------

# All datasets
# model = "GammaIBP_eb"; # var_GammaIBP = 1
model = "classicIBP_eb"
res_vec = apply(data, 1, function(counts){
  fit <- tryCatch(ParEst_PFFA(n,counts,model,var_GammaIBP  = var_GammaIBP ),
                  error=function(e) NA)
  fit
})


res_mat = matrix(0,nrow = BB, ncol = 0)
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$sigma_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$c_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$gam_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$fval) )
colnames(res_mat) <- names(res_vec[[1]])

colMeans(res_mat)

fitcheck_all = apply(res_mat, 1, function(x) {
  fk_stable_beta_process( c = x[2], sigma = x[1], 
                          gamma = x[3], base_rnd = base_rnd,
                          eps = 1e-6, seed = 12323123)$w 
})
fitcheck = lapply(fitcheck_all, function(x) x[1:xmax])
fitcheck = do.call(cbind, fitcheck); fitcheck[is.na(fitcheck)] = 0
fitcheck_qnt = apply(fitcheck, 1, quantile, probs = c(0.025,0.5,0.975))

## Plot data
ymax = max(max(qnt_3IBP),max(fitcheck_qnt)) * 1.02
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1)
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_3IBP[1,1:xmax], rev(qnt_3IBP[3,1:xmax])),
         col = "grey85",
         border = NA) # plot in-sample bands
polygon( c(1:xmax, rev(1:xmax)),
         c(fitcheck_qnt[1,1:xmax], rev(fitcheck_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands


# Sim. data from Zipfs -------------------------------------------------------
BB = 100
s = 0.95

M = 120
ptrue = sim_TruncatedZipfs_features(M = M, s = s)
# Gen. all datasets
data = lapply(ptrue, function(pj) rbinom(n = BB, size = n, prob = pj) )
data = do.call(cbind, data)

par(mfrow = c(1,1), mgp=c(1.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,type = "n", xlim = c(0,M), ylim = c(0,1), 
     xlab = "", ylab = "Prob.")
points( ptrue, pch = 16, cex = 0.5 )

## Fit params. Ale --------------------------------------------------------------
# All datasets
start_params <- c(alpha = 0.1, gamma= 1, u = 1)
res_vec = apply(data, 1, function(counts){
  counts = counts[which(counts > 0)]
  Kn = length(counts)
  fit <- tryCatch(optim(par = start_params, fn = llik_PP3Parm, 
                        method = "L-BFGS-B",
                        n = n, Kn = Kn, data_obs = counts,
                        lower = c(1e-16, 1e-16, 1e-16), 
                        upper = c(1-1e-10, Inf, Inf)) ,
                  error=function(e) NA)
  list("sigma_mle" = fit$par[1],
       "c_mle" = fit$par[3] - fit$par[1],
       "gamma_mle" = fit$par[2],
       "fval" = fit$value)
})


res_mat = matrix(0,nrow = BB, ncol = 0)
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$sigma_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$c_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$gamma_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$fval) )
colnames(res_mat) <- names(res_vec[[1]])

xmax = 100
fitcheck_all = apply(res_mat, 1, function(x) {
  fk_stable_beta_process( c = x[2], sigma = x[1], 
                          gamma = x[3], base_rnd = base_rnd,
                          eps = 1e-6, seed = 12323123)$w 
})
fitcheck = lapply(fitcheck_all, function(x) x[1:xmax])
fitcheck = do.call(cbind, fitcheck); fitcheck[is.na(fitcheck)] = 0
fitcheck_qnt = apply(fitcheck, 1, quantile, probs = c(0.025,0.5,0.975))

## Plot data
ymax = max(max(fitcheck_qnt)) * 1.02
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1)
polygon( c(1:xmax, rev(1:xmax)),
         c(fitcheck_qnt[1,1:xmax], rev(fitcheck_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points( ptrue, pch = 16, cex = 0.5 )

## Fit params --------------------------------------------------------------

# All datasets
# model = "GammaIBP_eb"; # var_GammaIBP = 1
model = "classicIBP_eb"
res_vec = apply(data, 1, function(counts){
  fit <- tryCatch(ParEst_PFFA(n,counts,model,var_GammaIBP  = var_GammaIBP ),
                  error=function(e) NA)
  fit
})


res_mat = matrix(0,nrow = BB, ncol = 0)
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$sigma_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$c_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$gam_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$fval) )
colnames(res_mat) <- names(res_vec[[1]])

colMeans(res_mat)

fitcheck_all = apply(res_mat, 1, function(x) {
  fk_stable_beta_process( c = x[2], sigma = x[1], 
                          gamma = x[3], base_rnd = base_rnd,
                          eps = 1e-6, seed = 12323123)$w 
})
fitcheck = lapply(fitcheck_all, function(x) x[1:xmax])
fitcheck = do.call(cbind, fitcheck); fitcheck[is.na(fitcheck)] = 0
fitcheck_qnt = apply(fitcheck, 1, quantile, probs = c(0.025,0.5,0.975))

## Plot data
ymax = max(max(fitcheck_qnt)) * 1.02
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1)
polygon( c(1:xmax, rev(1:xmax)),
         c(fitcheck_qnt[1,1:xmax], rev(fitcheck_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points( ptrue, pch = 16, cex = 0.5 )



