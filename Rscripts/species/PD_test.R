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
BB = 100

sigma_true = 0.1; theta_true = 50
xmax = 500

sim_PD_mat = r_SB(BB,xmax, sigma_true, theta_true,seed) # BB x xmax
sim_PD_mat = apply(sim_PD_mat, 1, sort, decreasing = TRUE) # xmax x BB
colSums(sim_PD_mat)
qnt_PD = apply(sim_PD_mat, 1, quantile, probs = c(0.025,0.5,0.975)) # 3 x xmax

## Plot data
ymax = max(qnt_PD) * 1.02
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,100), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_PD[1,1:xmax], rev(qnt_PD[3,1:xmax])),
         col = "grey85",
         border = NA) # plot in-sample bands

# Gen. all datasets
counts <- vapply( seq_len(BB), function(b) as.integer(rmultinom(1, size = n, prob = sim_PD_mat[, b])),
                  integer(xmax) ) # xmax x BB
dim(counts)
colSums(counts) 

## Fit and check Ale -----------------------------------------------------------
# All datasets
start_params <- c(alpha = 0.5, theta = 1)
res_vec = apply(counts, 2, function(Nj){
  Nj = Nj[which(Nj > 0)]
  Kn = length(Nj)
  fit <- tryCatch(optim(par = start_params, fn = llik_pyp, 
                        n = n, Kn = Kn, data_obs = Nj, # extra parameters
                        method = "L-BFGS-B",
                        lower = c(0, -1), upper = c(1-1e-10, Inf)) ,
                  error=function(e) NA)
  list("sigma_mle" = fit$par[1],
       "theta_mle" = fit$par[2],
       "fval" = fit$value)
})

res_mat = matrix(0,nrow = BB, ncol = 0)
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$sigma_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$theta_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$fval) )
colnames(res_mat) <- names(res_vec[[1]])

round(colMeans(res_mat),2)

fitcheck_all = apply(res_mat, 1, function(x) {
  r_SB(1,xmax, x[1], x[2], seed)
})
fitcheck_all = apply(fitcheck_all, 2, sort, decreasing = TRUE) 
fitcheck_qnt = apply(fitcheck_all, 1, quantile, probs = c(0.025,0.5,0.975))

## Plot data
ymax = max(max(qnt_PD),max(fitcheck_qnt)) * 1.02
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,100), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1)
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_PD[1,1:xmax], rev(qnt_PD[3,1:xmax])),
         col = "grey85",
         border = NA) # plot in-sample bands
polygon( c(1:xmax, rev(1:xmax)),
         c(fitcheck_qnt[1,1:xmax], rev(fitcheck_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands



