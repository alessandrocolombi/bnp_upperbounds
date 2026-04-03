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


# Plot options ------------------------------------------------------------------
save_img = FALSE
width = 12; height = 6
cex.labels <- cex.lab <- 2
cex.axis <- 2
cex.legend <- 1.5
mycol = c("darkgreen","darkorange","darkred","darkblue","lightblue")
mycol2 = c("black","lightblue")
lgd_names = c("Unbounded","Bounded","IBP","MBP","FB")
xmax = 140
# Load --------------------------------------------------------------------

load("RawDataInc.Rdat")
n = ncol(A)
Kn = nrow(A)
N_j = rowSums(A)
names(N_j) = as.character(1:Kn)
Nj_ordered = sort(N_j, decreasing = TRUE)

seed = 34231
set.seed(seed)

data = t(A)


# Upper bounds ------------------------------------------------------------
alfa <- alpha <- 0.05; 
var_fct = 500
Rmax = 100

## Freq. ------------------------------------------------------------
# Bounded
b_n <- log(n)
Mguess = 50 + Kn
Nj_guess = c(N_j, rep(0,Mguess - length(N_j) ))
ubFreqBdd <- compute_UB_analytical(n, Nj_guess, Mguess, b_n, alfa, FALSE)
ubFreqBdd
# Unbounded
beta = 1e-5
Shat  <- sum(N_j) / n
Sstar <- ( sqrt( -log(beta) / (2 * n) ) +
             sqrt( Shat + (-log(beta) / (2 * n)) ) )^2
r_n   <- log( Sstar / (-log(1 - alfa + beta)) ) + log(n) - log(log(n))
ubFreqUbd <- compute_UB_rnorm(n, alfa, beta, r_n, Shat)
ubFreqUbd
## 3IBP ------------------------------------------------------------
# Param. estimation (3 params PP)
start_params <- c(alpha = 0.1, gamma= 1, u = 1)
fit <- optim(par = start_params, fn = llik_PP3Parm, 
             method = "L-BFGS-B",
             n = n, Kn = Kn, data_obs = N_j,
             lower = c(1e-16, 1e-16, 1e-16), 
             upper = c(1-1e-10, Inf, Inf)) 
alpha_mle = fit$par[1]
gamma_mle = fit$par[2]
c_mle     = fit$par[3] - alpha_mle

Nrep = 100
seeds = sample(1:999999, size = Nrep)
sim_3IBP = lapply(seeds, function(seed) fk_stable_beta_process( c = c_mle, sigma = alpha_mle, gamma = gamma_mle,
                                                                base_rnd = base_rnd, eps = 1e-6, seed = seed)$w )
sim_3IBP = lapply(sim_3IBP, function(x) x[1:xmax])
mat_3IBP = do.call(cbind, sim_3IBP)
qnt_3IBP = apply(mat_3IBP, 1, quantile, probs = c(0.025,0.5,0.975))

ymax = max(max(Nj_ordered)/n , max(qnt_3IBP) )
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

if(save_img)
  pdf("img/Crim4features_IBP_paramfit.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(2.25,0.5,0), mar = c(2,3.25,1,1), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "Prob. (IBP)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_3IBP[1,1:xmax], rev(qnt_3IBP[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:Kn, y = Nj_ordered/n, col = "darkred", pch = 16, cex = 0.5)
# segments(x0 = (1:Kn), x1 = (1:Kn), 
# y0 = rep(0,Kn), y1 = Nj_ordered/n, col = "darkred", lwd = 3)
points(x = (1:xmax)+0.25, y = qnt_3IBP[2,1:xmax], 
       col = "black", pch = 4)
# segments(x0 = (1:xmax)+0.25, x1 = (1:xmax)+0.25,
# y0 = rep(0,xmax), y1 = qnt_3IBP[2,1:xmax], 
# col = "black", lwd = 1, lty = 2)
if(save_img)
  dev.off()


# Upper bound (3 params PP)
ub3IBP = exp(compute_log_UBMarkov_BeBePois(Rmax, alpha_mle, c_mle, gamma_mle, n, alfa ))
ub3IBP

## MBP ---------------------------------------------------------------------
eb_init_BB <- list(alpha = -1, s = 100, Nhat_prime = 50)
eb_known_BB <- list()
eb_params_obj_BB <- eb_params(model = "BB", init = eb_init_BB, known = eb_known_BB )

res = GibbsFA_eb(feature_matrix = data,
                 model = "NegBinBB_eb", type = "EFPF",
                 eb_params =  eb_params_obj_BB, 
                 var_fct = var_fct)

a_mle = res$alpha+1
b_mle = res$theta - a_mle
r_nb = res$n0
q_nb = 1 - 1/res$var_fct; p_nb = 1/var_fct
kappa_n = exp( lgamma(b_mle+n) + lgamma(a_mle+b_mle) - lgamma(b_mle) - lgamma(a_mle+b_mle+n) )


M_mc = rnbinom(n = 10000, size = r_nb, prob = p_nb)
w_mc_list = lapply(M_mc, function(m) {w = rbeta(n = m, shape1 = a_mle, shape2 = b_mle); w = sort(w,decreasing = TRUE); w[1:xmax]})
mat_MixBin = do.call(cbind, w_mc_list)
mat_MixBin[is.na(mat_MixBin)] = 0
qnt_MixBin = apply(mat_MixBin, 1, quantile, probs = c(0.025,0.5,0.975), na.rm = TRUE)

ymax = max(max(Nj_ordered)/n , max(qnt_MixBin) )
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

if(save_img)
  pdf("img/Crim4features_MBP_paramfit.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(2.25,0.5,0), mar = c(2,3.25,1,1), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "Prob. (MBP)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_MixBin[1,1:xmax], rev(qnt_MixBin[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:Kn, y = Nj_ordered/n, col = "darkred", pch = 16)
# segments(x0 = (1:Kn), x1 = (1:Kn), 
#          y0 = rep(0,Kn), y1 = Nj_ordered/n, col = "darkred", lwd = 3)
points(x = (1:xmax)+0.25, y = qnt_MixBin[2,1:xmax], 
       col = "black", pch = 4)
# segments(x0 = (1:xmax)+0.25, x1 = (1:xmax)+0.25,
#          y0 = rep(0,xmax), y1 = qnt_MixBin[2,1:xmax], 
#          col = "black", lwd = 1, lty = 2)
if(save_img)
  dev.off()


ubMixBin = exp(compute_log_UBMarkov_BeBeMixNBin( Rmax, a_mle, b_mle, n, Kn, r_nb, p_nb, alpha))
ubMixBin = min(ubMixBin,1)

## FB ----------------------------------------------------------------------
Mguess = 50 + Kn
Nj_guess = c(N_j, rep(0,Mguess - length(N_j) ))
start_params <- c(a = 1, b = 1)
fit <- optim(par = start_params, fn = llik_FB,
             method = "L-BFGS-B",
             n = n, Kn = Kn, data_obs = Nj_guess, M=Mguess,
             lower = c(1e-10, 1e-10), upper = c(Inf, Inf))
a_FB = fit$par[1]; 
b_FB = fit$par[2]; 

w_mc_list = lapply(rep(Mguess,1000), function(m) {w = rbeta(n = m, shape1 = a_mle, shape2 = b_mle); w = sort(w,decreasing = TRUE); w[1:xmax]})
mat_FB = do.call(cbind, w_mc_list)
mat_FB[is.na(mat_FB)] = 0
qnt_FB = apply(mat_FB, 1, quantile, probs = c(0.025,0.5,0.975), na.rm = TRUE)

ymax = max(max(Nj_ordered)/n , max(qnt_FB) )
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

if(save_img)
  pdf("img/Crim4features_FB_paramfit.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(2.25,0.5,0), mar = c(2,3.25,1,1), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "Probs. (FB)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_FB[1,1:xmax], rev(qnt_FB[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:Kn, y = Nj_ordered/n, col = "darkred", pch = 16)
segments(x0 = (1:Kn), x1 = (1:Kn), 
         y0 = rep(0,Kn), y1 = Nj_ordered/n, col = "darkred", lwd = 3)
points(x = (1:xmax)+0.25, y = qnt_FB[2,1:xmax], 
       col = "black", pch = 4)
segments(x0 = (1:xmax)+0.25, x1 = (1:xmax)+0.25,
         y0 = rep(0,xmax), y1 = qnt_FB[2,1:xmax], 
         col = "black", lwd = 1, lty = 2)
if(save_img)
  dev.off()

## ii) Upper bound computation (FB)
ubFB = exp(compute_log_UBMarkov_FB( Rmax, a_FB, b_FB, n, Kn, Mguess, alpha))
ubFB = min(ubFB,1)

# Missing mass analysis ------------------------------------------------------------------
## Freq ------------------------------------------------------------
Chao = SpadeR::ChaoSpecies( A, datatype = "incidence_raw")
Chao
## Mixed-Binomial ------------------------------------------------------------
post_size = r_nb + Kn
post_Prfail = q_nb*kappa_n
ExpMstar = post_size * (post_Prfail)/(1 - post_Prfail)
ExpMstar

Mstar_mc = rnbinom(n = 10000,size = post_size, prob = 1 - post_Prfail)

pMstar = rep(0,max(Mstar_mc))
for(k in 1:(length(pMstar)-1)){
  pMstar[k] = dnbinom(x = (k-1), size = post_size, prob = 1 - post_Prfail)
}

bp1 = barplot(pMstar)

ypos = seq(0, max(pMstar), length.out = 3)
ynames = round(ypos,2)


if(save_img)
  pdf("img/Crim4features_MixBin_Mstar.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(2,0.5,0), mar = c(1,3.5,1,0), cex = 1.5)
barplot( height = pMstar, 
         names.arg = "", las = 2, col = "darkred", border = NA,
         main = " ", ylab = expression("Probs. " * (M^"*" )), yaxt = "n" )
axis( side = 2, at = ypos, labels = ynames ,las = 1)
text( x = bp1, y = -0.001, 
      labels = 0:(length(pMstar)-1), 
      srt = 0, adj = 0.5, xpd = TRUE, cex = 1)
if(save_img)
  dev.off()

