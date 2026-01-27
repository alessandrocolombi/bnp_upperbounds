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
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")

# Plot options ------------------------------------------------------------------
save_img = FALSE
width = 12; height = 6
cex.labels = 2
xmax = 52
# Load --------------------------------------------------------------------

load("RawDataInc.Rdat")
# ni = 50
# A = A[sample(1:n, size = ni),]
n = nrow(A)
Kn = ncol(A)
N_j = colSums(A)
names(N_j) = as.character(1:Kn)
Nj_ordered = sort(N_j, decreasing = TRUE)
seed = 34231
set.seed(seed)

# Upper bounds ------------------------------------------------------------
alfa = 0.05; 
Rmax = 100
var_gamma = 10
var_nb    = 1000

## Freq. ------------------------------------------------------------
# Bounded
b_n <- log(n)
Mguess = 100 + Kn
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
  pdf("img/Mod2_Inc4meet_3IBP_paramfit.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_3IBP[1,1:xmax], rev(qnt_3IBP[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:Kn, y = Nj_ordered/n, col = "darkred", pch = 16)
segments(x0 = (1:Kn), x1 = (1:Kn), 
         y0 = rep(0,Kn), y1 = Nj_ordered/n, col = "darkred", lwd = 3)
points(x = (1:xmax)+0.25, y = qnt_3IBP[2,1:xmax], 
       col = "black", pch = 4)
segments(x0 = (1:xmax)+0.25, x1 = (1:xmax)+0.25,
         y0 = rep(0,xmax), y1 = qnt_3IBP[2,1:xmax], 
         col = "black", lwd = 1, lty = 2)
if(save_img)
  dev.off()

# Upper bound (3 params PP)
ub3IBP = exp(compute_log_UBMarkov_BeBePois(Rmax, alpha_mle, c_mle, gamma_mle, n, alfa ))
ub3IBP

## Mixed-Poisson ------------------------------------------------------------

# Param. estimation (Mixed Poisson)
start_params <- c(alpha = 0.1, u = 1, mu_gamma = 1)
fit <- optim(par = start_params, fn = llik_MixPois,
             method = "L-BFGS-B",
             n = n, Kn = Kn, data_obs = N_j, var_gamma = var_gamma,
             lower = c(1e-16, 1e-16, 1e-16), upper = c(1-1e-10, Inf, Inf))
alpha_mle = fit$par[1]
c_mle = fit$par[2] - alpha_mle
mugamma_mle = fit$par[3]

gamma_hyperparams = gamma_shape_rate(mugamma_mle,var_gamma)
u = gamma_hyperparams$shape
v = gamma_hyperparams$rate
# Upper bound (Mixed Poisson)
ubMixPois =  exp(compute_log_UBMarkov_BeBeMixPois( Rmax, alpha_mle, c_mle, n, Kn, u, v, alfa))
ubMixPois

## Mixed Binomial ------------------------------------------------------------
# Param. estimation (Mixed Binomial)
start_params <- c(a = 1, b = 1, mu_nb = 10)
fit <- optim(par = start_params, fn = llik_MixBin,
             method = "L-BFGS-B",
             n = n, Kn = Kn, data_obs = N_j, var_nb = var_nb,
             lower = c(1e-10, 1e-10, 1e-10), upper = c(Inf, Inf, var_nb-1e-10))
a_mle = fit$par[1]
b_mle = fit$par[2] 
munb_mle = fit$par[3]

nb_hyperparams = NegBin_params(munb_mle,var_nb)
r_nb = nb_hyperparams$r
p_nb = nb_hyperparams$p
q_nb = 1 - p_nb
# Upper bound (Mixed Binomial)
ubMixBin = exp(compute_log_UBMarkov_BeBeMixNBin( Rmax, a_mle, b_mle, n, Kn, r_nb, p_nb, alfa))
ubMixBin



M_mc = rnbinom(n = 50000, size = r_nb, prob = p_nb)
w_mc_list = lapply(M_mc, function(m) {w = rbeta(n = m, shape1 = a_mle, shape2 = b_mle); w = sort(w,decreasing = TRUE); w[1:xmax]})
mat_MixBin = do.call(cbind, w_mc_list)
mat_MixBin[is.na(mat_MixBin)] = 0
qnt_MixBin = apply(mat_MixBin, 1, quantile, probs = c(0.025,0.5,0.975), na.rm = TRUE)

ymax = max(max(Nj_ordered)/n , max(qnt_MixBin) )
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

if(save_img)
  pdf("img/Mod2_Inc4meet_MixBin_paramfit.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_MixBin[1,1:xmax], rev(qnt_MixBin[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:Kn, y = Nj_ordered/n, col = "darkred", pch = 16)
segments(x0 = (1:Kn), x1 = (1:Kn), 
         y0 = rep(0,Kn), y1 = Nj_ordered/n, col = "darkred", lwd = 3)
points(x = (1:xmax)+0.25, y = qnt_MixBin[2,1:xmax], 
       col = "black", pch = 4)
segments(x0 = (1:xmax)+0.25, x1 = (1:xmax)+0.25,
         y0 = rep(0,xmax), y1 = qnt_MixBin[2,1:xmax], 
         col = "black", lwd = 1, lty = 2)
if(save_img)
  dev.off()



## Summary -----------------------------------------------------------------
c(ub3IBP,ubMixPois,ubMixBin,ubFreqBdd,ubFreqUbd)


# Missing mass analysis ------------------------------------------------------------------

## Freq ------------------------------------------------------------
# Chao = SpadeR::ChaoSpecies(t(A),datatype = "incidence_raw")
# Fails because there are no singleton or doubletons
## 3IBP ------------------------------------------------------------
## Mixed-Poisson ------------------------------------------------------------
## Mixed-Binomial ------------------------------------------------------------
kappa_n = exp( lgamma(b_mle+n) + lgamma(a_mle+b_mle) - lgamma(b_mle) - lgamma(a_mle+b_mle+n) )
post_size = r_nb + Kn
post_Prfail = q_nb*kappa_n
ExpMstar = post_size * (post_Prfail)/(1 - post_Prfail)
ExpMstar

Mstar_mc = rnbinom(n = 10000,size = post_size, prob = 1 - post_Prfail)
table(Mstar_mc)/10000

pMstar = rep(0,5)
for(k in 1:(length(pMstar)-1)){
  pMstar[k] = dnbinom(x = (k-1), size = post_size, prob = 1 - post_Prfail)
}

bp1 = barplot(pMstar)

if(save_img)
  pdf("img/Mod2_Inc4meet_MixBin_Mstar.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(2,0.5,0), mar = c(2,3,1,0), cex = 2)
barplot( height = pMstar, 
         names.arg = "", las = 2, col = "darkred", border = NA,
         main = " ", ylab = "Mstar", yaxt = "n" )
axis( side = 2, at = seq(0, max(pMstar), by=0.1), las = 1)
text( x = bp1, y = -0.1, 
      labels = 0:(length(pMstar)-1), 
      srt = 0, adj = 0.5, xpd = TRUE)
if(save_img)
  dev.off()









