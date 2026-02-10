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

# Common params -----------------------------------------------------------
seed = 42
set.seed(seed)
alfa <- alpha <- 0.05
n = 100

# Sim. data from model -------------------------------------------------------
BB = 10

mu = 120;v  = 240
r_nb = (mu*mu)/(v-mu)
p_nb = mu/v
q_nb = 1-p_nb
a_beta = 0.5; b_beta = 10
a_beta/(a_beta+b_beta)

M_mc = rnbinom(n = BB, size = r_nb, prob = p_nb)
MM = max(M_mc)

ptrue = matrix(0,nrow = BB, ncol = MM)
for(b in 1:BB){
  ptrue[b,1:M_mc[b]] = rbeta(n=M_mc[b],shape1 = a_beta, shape2 = b_beta)
}

ptrue_sorted = apply(ptrue, 1, sort, decreasing = TRUE)
qnt_Ptrue = apply(ptrue_sorted, 1, quantile, probs = c(0.025,0.5,0.975), na.rm = TRUE)

# Gen. all datasets
data = apply(ptrue, 1 ,function(pj) rbinom(n = MM, size = n, prob = pj) )
data = t(data)
dim(data)


## Fit and check -----------------------------------------------------------

# Single dataset
counts = data[1,]
model = "NegBinBB_eb"
var_fct_vec = c(1,5,10,50,100,500,1000)
# var_fct_vec = c(100)
# var_fct<-x <- 100
res_vec = lapply(data, function(x){
  fit <- tryCatch(ParEst_PFFA(n,counts,model,var_fct = x),
                  error=function(e) NA)
  fit
})

sapply(res_vec, function(x) x$fval)
sapply(res_vec, function(x) x$r_nb_mle)
sapply(res_vec, function(x) x$mu_nb_mle)
sapply(res_vec, function(x) x$a_mle)
sapply(res_vec, function(x) x$b_mle)

# All datasets
var_fct <- 100
res_vec = apply(data, 1, function(counts){
  fit <- tryCatch(ParEst_PFFA(n,counts,model,var_fct = var_fct),
                  error=function(e) NA)
  fit
})

res_mat = matrix(0,nrow = BB, ncol = 0)
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$a_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$b_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$r_nb_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$mu_nb_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$p_nb_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$fval) )
colnames(res_mat) <- names(res_vec[[1]])

w_mc_list = apply(res_mat, 1, function(x) {
  m = rnbinom(n = 1, size = x[3], prob = x[5])
  w = rbeta(n = m, shape1 = x[1], shape2 = x[2]) 
  w = sort(w,decreasing = TRUE) 
  w
})
Mmax = median(sapply(w_mc_list, length))
w_mc_list = lapply(w_mc_list, function(x) {
  y = rep(0,Mmax); m = length(x);
  y[1:min(Mmax,m)] = x[1:min(Mmax,m)]; y}
)
mat_MixBin = do.call(cbind, w_mc_list)
mat_MixBin[is.na(mat_MixBin)] = 0
qnt_MixBin = apply(mat_MixBin, 1, quantile, probs = c(0.025,0.5,0.975), na.rm = TRUE)


ymax = max(max(ptrue), max(qnt_MixBin) )
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,Mmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:Mmax, rev(1:Mmax)),
         c(qnt_MixBin[1,1:Mmax], rev(qnt_MixBin[3,1:Mmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
polygon( c(1:MM, rev(1:MM)),
         c(qnt_Ptrue[1,1:MM], rev(qnt_Ptrue[3,1:MM])),
         col = "grey85",
         border = NA) # plot in-sample bands

# Sim. data from Zipfs -------------------------------------------------------
BB = 100
s = 1.25

M = 120
ptrue = sim_TruncatedZipfs_features(M = M, s = s)
# Gen. all datasets
data = lapply(ptrue, function(pj) rbinom(n = BB, size = n, prob = pj) )
data = do.call(cbind, data)

par(mfrow = c(1,1), mgp=c(1.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,type = "n", xlim = c(0,M), ylim = c(0,1), 
     xlab = "", ylab = "Prob.")
points( ptrue, pch = 16, cex = 0.5 )

## Fit params --------------------------------------------------------------

# Single dataset
counts = data[1,]
model = "NegBinBB_eb"
var_fct_vec = c(1,5,10,50,100,500,1000)
# var_fct_vec = c(100)
# var_fct<-x <- 100
res_vec = lapply(data, function(x){
  fit <- tryCatch(ParEst_PFFA(n,counts,model,var_fct = x),
                  error=function(e) NA)
  fit
})

sapply(res_vec, function(x) x$fval)
sapply(res_vec, function(x) x$r_nb_mle)
sapply(res_vec, function(x) x$mu_nb_mle)
sapply(res_vec, function(x) x$a_mle)
sapply(res_vec, function(x) x$b_mle)

# All datasets
var_fct <- 100
res_vec = apply(data, 1, function(counts){
  fit <- tryCatch(ParEst_PFFA(n,counts,model,var_fct = var_fct),
                  error=function(e) NA)
  fit
})

res_mat = matrix(0,nrow = BB, ncol = 0)
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$a_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$b_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$r_nb_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$mu_nb_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$p_nb_mle) )
res_mat = cbind(res_mat, sapply(res_vec, function(x) x$fval) )
colnames(res_mat) <- names(res_vec[[1]])

w_mc_list = apply(res_mat, 1, function(x) {
  m = rnbinom(n = 1, size = x[3], prob = x[5])
  w = rbeta(n = m, shape1 = x[1], shape2 = x[2]) 
  w = sort(w,decreasing = TRUE) 
  w
})
Mmax = 50;#median(sapply(w_mc_list, length))
w_mc_list = lapply(w_mc_list, function(x) {
  y = rep(0,Mmax); m = length(x);
  y[1:min(Mmax,m)] = x[1:min(Mmax,m)]; y}
)
mat_MixBin = do.call(cbind, w_mc_list)
mat_MixBin[is.na(mat_MixBin)] = 0
qnt_MixBin = apply(mat_MixBin, 1, quantile, probs = c(0.025,0.5,0.975), na.rm = TRUE)


ymax = max(max(ptrue), max(qnt_MixBin) )
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,Mmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:Mmax, rev(1:Mmax)),
         c(qnt_MixBin[1,1:Mmax], rev(qnt_MixBin[3,1:Mmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points( ptrue, pch = 16, cex = 0.5 )

