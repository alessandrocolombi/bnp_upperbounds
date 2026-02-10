# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/species/")
setwd(wd)

# Functions ---------------------------------------------------------------
source("../../R/Rfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")

# Funzioni ----------------------------------------------------------------

# Common params -----------------------------------------------------------
seed = 42
set.seed(seed)
alfa <- alpha <- 0.05
n = 1000


# Sim. data from model -------------------------------------------------------
BB = 50
M = 100
gamma_true = 0.15
ptrue = matrix(0,nrow = M, ncol = BB)
ptrue = apply(ptrue, 2, function(x) {w = rgamma(n = M, shape = gamma_true, rate = 1); w/sum(w)})
ptrue = apply(ptrue, 2, sort, decreasing = TRUE) 
colSums(ptrue)
qnt_FDP = apply(ptrue, 1, quantile, probs = c(0.025,0.5,0.975)) # 3 x xmax

## Plot data
ymax = max(qnt_FDP) * 1.02
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,M), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:M, rev(1:M)),
         c(qnt_FDP[1,1:M], rev(qnt_FDP[3,1:M])),
         col = "grey85",
         border = NA) # plot in-sample bands

# Gen. all datasets
counts <- vapply( seq_len(BB), function(b) as.integer(rmultinom(1, size = n, prob = ptrue[, b])),
                  integer(M) ) 
dim(counts)
colSums(counts) 

## Fit and check Ale -----------------------------------------------------------
# All datasets
start_params <- c(gamma = 0.5)
res_vec = apply(counts, 2, function(Nj){
  fit <- tryCatch(optim(par = start_params, fn = llik_DirMult, 
                        n = n, M = M, data = Nj, # extra parameters
                        method = "L-BFGS-B",
                        lower = c(1e-10), upper = c(Inf)) ,
                  error=function(e) NA)
  list("gamma_mle" = fit$par[1],
       "fval" = fit$value)
})

res_mat = matrix(0,nrow = BB, ncol = 0)
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$gamma_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$fval) )
colnames(res_mat) <- names(res_vec[[1]])

round(colMeans(res_mat),2)

fitcheck_all = apply(res_mat, 1, function(x) {
  w = rgamma(n = M, shape = x[1], rate = 1); w/sum(w)
})
fitcheck_all = apply(fitcheck_all, 2, sort, decreasing = TRUE) 
fitcheck_qnt = apply(fitcheck_all, 1, quantile, probs = c(0.025,0.5,0.975))

## Plot data
ymax = max(max(qnt_FDP),max(fitcheck_qnt)) * 1.02
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,M), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1)
polygon( c(1:M, rev(1:M)),
         c(qnt_FDP[1,1:M], rev(qnt_FDP[3,1:M])),
         col = "grey85",
         border = NA) # plot in-sample bands
polygon( c(1:M, rev(1:M)),
         c(fitcheck_qnt[1,1:M], rev(fitcheck_qnt[3,1:M])),
         col = "#FDE333",
         border = NA) # plot in-sample bands


# Upper bound -------------------------------------------------------------
Nj = counts[,1]
Kn = length(which(Nj > 0))
Rmax = 100
# Upper bound (DirMulti)
ubDM = exp(compute_log_UB_DirMulti( Rmax, gamma_true, M, Kn, n, alpha ))
ubDM = min(ubDM,1)
ubDM

