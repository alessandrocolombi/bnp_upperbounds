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

# Colors ------------------------------------------------------------------


# Load --------------------------------------------------------------------

load("RawDataInc.Rdat")
n = nrow(A)
Kn = ncol(A)
N_j = colSums(A)
names(N_j) = as.character(1:Kn)

seed = 34231
set.seed(seed)

# Upper bounds ------------------------------------------------------------
alfa = 0.05; 
Rmax = 100
var_gamma = 10
var_nb    = 10

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
start_params <- c(a = 1, b = 1, mu_nb = 1)
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
# Upper bound (Mixed Binomial)
ubMixBin = exp(compute_log_UBMarkov_BeBeMixNBin( Rmax, a_mle, b_mle, n, Kn, r_nb, p_nb, alfa))
ubMixBin


# Missing mass analysis ------------------------------------------------------------------

## Freq ------------------------------------------------------------
# Fails because there are no singleton or doubletons
## 3IBP ------------------------------------------------------------
## Mixed-Poisson ------------------------------------------------------------
## Mixed-Binomial ------------------------------------------------------------
