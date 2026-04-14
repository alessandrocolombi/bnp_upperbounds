# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc, wd_unicatt, wd_g100, wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd, "bnp_upperbounds/Rscripts/species")
setwd(wd)

compute_ab_beta = function(m, v, tol = sqrt(.Machine$double.eps) ) 
{
  if (length(m) != 1 || length(v) != 1) 
    stop("m and v must be scalars")
  if (!is.finite(m) || !is.finite(v)) 
    stop("m and v must be finite")
  if (m <= 0 || m >= 1) 
    stop("m must satisfy 0 < m < 1")
  if (v <= 0) 
    stop("v must satisfy v > 0")
  
  vmax = m * (1 - m)
  
  if (v >= vmax) 
    stop("Need v < m*(1-m) for a proper Beta distribution")
  
  # detect near-boundary regime: a+b very close to 0
  if ((vmax - v) <= tol * vmax) {
    warning("v is extremely close to m*(1-m): a and b are near 0 and numerically unstable")
  }
  
  kappa = (vmax - v) / v
  a = m * kappa
  b = (1 - m) * kappa
  
  c(a = a, b = b)
}
lpost_FB <- function(x, n, Kn, M, data, hy){
  mu <- x[1]; kappa <- x[2]
  var <- (mu*(1-mu))/kappa
  if( mu < 1e-8 || mu > (1-1e-8) ){
    return( +100000000 )
  }
  if( kappa < 1e-8 ){
    return( +100000000 )
  }
  a_mu  <- hy[1]; b_mu <- hy[2]; a_kappa <- hy[3]; b_kappa <- hy[4] 
  log_prior_mu     = (a_mu - 1)*log(mu) + (b_mu - 1)*log(1-mu)
  log_prior_kappa  = (a_kappa - 1)*log(kappa)   - b_kappa*kappa
  ab = compute_ab_beta(mu,var)
  llik = log_efpf_FB(n,Kn,M,data,ab[1],ab[2])
  return( -(llik+log_prior_mu+log_prior_kappa) )
  
} 



## FB hyperparameters
a_mu = 1; b_mu = 1
mu_kappa = 100; var_kappa = 10
a_kappa = mu_kappa*mu_kappa/var_kappa
b_kappa = mu_kappa/var_kappa
hy_FB = c(a_mu,b_mu,a_kappa,b_kappa)

# fit
start_params <- c(mu = 0.1, kappa = 100)
fit <- optim(par = start_params, fn = lpost_FB,
             method = "L-BFGS-B",
             n=n, Kn=Kn, data=n_i, M=M, hy=hy_FB,
             lower = c(1e-6, 1e-6), upper = c(1-1e-6, Inf))
mu_est    = fit$par[1]
kappa_est = fit$par[2]
var_est <- (mu_est*(1-mu_est))/kappa_est
ab_est = compute_ab_beta(mu_est,var_est)
ab_est