
# Gen. Dist. (species) ----------------------------------------------------
sim_zipfs = function(M,s){
  w = sapply(1:M,function(j) j^(-s))
  w / sum(w)
}
sim_geom = function(M,a ){
  w = sapply(0:(M-1),function(j) (1-a)^(j-1)  )
  w / sum(w)
}
sim_negbin = function(M,l,r){
  w = dnbinom(x = 0:(M-1), size = l, prob = r) # Painsky parametrization
  w / sum(w)
}
sim_betabin = function(M,a,b){
  w = dbetabinom.ab(x = 0:(M-1), size = M-1, shape1 = a, shape2 = b)
  # lw = sapply(1:M,function(j) lchoose(M,j) + lbeta(j+a,M-j+b) - lbeta(a,b) )
  # lw = lchoose(M,1:M) + lbeta(1:M + a, M - 1:M + b )  
  # max_lw = max(lw)
  # w = exp(lw - max_lw)*exp(max_lw)
  w / sum(w)
}
sim_unif = function(M){
  w = 1:M
  w / sum(w)
}

# Upper bound Painsky -----------------------------------------------------
ub_pain = function(n,Rmax,alfa){
  rgrid = 1:Rmax
  qstar_pan = function(n,r){
    (r-1)/(r-1+n)
  }
  pan_ub_grid = sapply(rgrid, function(r){
    qstar = qstar_pan(n,r)
    res = qstar^(r-1) * (1 - qstar )^n * (1/alfa)
    res^(1/r)
  })
  pan_ub = min(pan_ub_grid)
  pan_ub
}


# Worst case dist. --------------------------------------------------------
compute_r_unbounded = function(n,Rmax,alfa){
  rgrid = 1:Rmax
  qstar_pan = function(n,r){
    (r-1)/(r-1+n)
  }
  pan_ub_grid = sapply(rgrid, function(r){
    qstar = qstar_pan(n,r)
    res = qstar^(r-1) * (1 - qstar )^n * (1/alfa)
    res^(1/r)
  })
  which.min(pan_ub_grid)
}


# Objective functions (species) -----------------------------------------------------
llik_pyp <- function(x, n, Kn, data_obs) {
  alpha <- x[1]
  theta <- x[2]
  -log_eppfPYP( n, Kn, data_obs, alpha, theta )
}
llik_FD <- function(x, n, Kn, data_obs, M_max) {
  gamma <- x[1]
  Lambda <- x[2]
  -log_eppfFD(n,Kn,data_obs,gamma,Lambda,M_max)
}

# Objective functions (features) -----------------------------------------------------
llik_PP <- function(x,n, Kn, data_obs, gamma) {
  alpha <- x[1]
  u <- x[2]
  c <- u - alpha
  -log_efpfBeBePois(n, Kn, data_obs, alpha, c, gamma )
}
llik_PP2 <- function(x, n, Kn, data_obs) {
  alpha <- 0 
  gamma <- x[1]
  c <- x[2]
  -log_efpfBeBePois(n, Kn, data_obs, alpha, c, gamma )
}
llik_PP3Parm <- function(x, n, Kn, data_obs) {
  alpha <- x[1]
  gamma <- x[2]
  u     <- x[3]
  c     <- u - alpha
  -log_efpfBeBePois(n, Kn, data_obs, alpha, c, gamma )
}
llik_MixPois <- function(x, n, Kn, data_obs, var_gamma) {
  alpha <- x[1]
  u <- x[2]
  c <- u - alpha
  mu_gamma <- x[3]
  -log_efpfBeBeMixPois( n, Kn, data_obs, alpha, c, mu_gamma, var_gamma )
}
llik_MixBin <- function(x, n, Kn, data_obs, var_nb) {
  a <- x[1]
  b <- x[2]
  mu_nb <- x[3]
  -log_efpfBeBeMixNBin( n, Kn, data_obs, a, b, mu_nb, var_nb )
}
f_beta <- function(x, n, alfa, Shat){
  if(n <= 0)
    stop("Error in lf_beta: n<=0")
  if(alfa-x <= 0)
    stop("Error in lf_beta: alfa-x<=0")
  if( 1/x <= 0 )
    stop("Error in lf_beta: 1/x<=0")
  n/(alfa-x) * (  sqrt( (log(1/x))/(2*n) ) + sqrt( (log(1/x))/(2*n) + Shat) )
}
lf_beta <- function(x, n, alfa, Shat){
  if(n <= 0)
    stop("Error in lf_beta: n<=0")
  if(alfa-x <= 0)
    stop("Error in lf_beta: alfa-x<=0")
  if( 1/x <= 0 )
    stop("Error in lf_beta: 1/x<=0")
  
  log(n) -log(alfa-x) + log(  sqrt( (log(1/x))/(2*n) ) + sqrt( (log(1/x))/(2*n) + Shat) )
}

# Gen. Dist. (features) ----------------------------------------------------
sim_zipfs_features = function(M,s){
  w = sapply(1:(M),function(j) j^(-s))
  w 
}
sim_geom_features = function(M,a){
  w = sapply(0:(M-1),function(j) (1-a)^(j-1)  )
}
sim_ghilo_features = function(M){
  m = floor(M/3)
  w = rep(0,M)
  w[1:m] = 0.015
  w[(m+1):(2*m)] = 0.01
  w[(2*m+1):M] = 0.005
  w 
}

# Mmax-based stopping rules -----------------------------------------------
SRabu_grid_multiple_run <- function(eps, data, nstart, seed0, Nrep, alpha, M_max)
{
  ## Functions
  suppressWarnings(suppressPackageStartupMessages(library(tibble)))
  suppressWarnings(suppressPackageStartupMessages(library(parallel)))
  suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
  suppressWarnings(suppressPackageStartupMessages(library(progress)))
  suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
  source("../../R/Rfunctions.R")
  Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
  
  
  set.seed(seed0)
  seeds = sample(1:999999, size = Nrep)
  res = matrix(nrow = Nrep, ncol = 3)
  colnames(res) = c("FD","PYP","Freq")
  
  for(ii in 1:Nrep){
    seed = seeds[ii]
    res[ii,] = SRabu_grid_single_run(eps, data, nstart, seed, alpha, M_max)
  }
  return(res)
}

SRabu_grid_single_run <- function(eps, data, nstart, seed, alpha, M_max)
{
  set.seed(seed) # set seed
  RmaxFD = 50; Rmax = 100
  nj = table(data) # compute abundances
  n  = sum(nj) # total num. obs.
  Kn = length(nj) # total num. distinct
  ordered_idx = sample(1:n, size = n) # choose ordering of obs.
  
  if(nstart >= (n-1))
    stop("nstart must be smaller than n-1")
  
  # Stopping flags and outputs
  stopped_FD   <- FALSE
  stopped_PYP  <- FALSE
  stopped_Freq <- FALSE
  Nstop_FD     <- NA_integer_
  Nstop_PYP    <- NA_integer_
  Nstop_Freq   <- NA_integer_
  
  ## ------------------------------------------------------------
  ## Run loop up to n_max = n
  ## ------------------------------------------------------------
  n_max = n
  ni = 2
  for(ni in nstart:(n_max-1)) {
    # Allow for a non-multiple n_max if needed
    remaining <- n_max - ni
    if (remaining <= 0L) break
    
    ## ---- Observed abundance vector (true + error species) ----
    idx_species_i = ordered_idx[1:ni] # select obs. up to time ni
    data_i = data[idx_species_i] 
    Nj_i = table(data_i) # compute frequencies
    Nj_i = Nj_i[Nj_i > 0]
    Kobs_i = length(which(Nj_i > 0))
    if( Kobs_i == 0L) next   # nothing observed yet
    
    ## ---- 5.1 FD (on observed data) ----
    if (!stopped_FD) {
      # Param. estimation (FD)
      start_params <- c(gamma = 0.1, Lambda = Kobs_i)
      fit <- optim(par = start_params, fn = llik_FD,
                   n = sum(Nj_i), Kn = Kobs_i, data_obs = Nj_i, M_max = M_max,# extra parameters
                   method = "L-BFGS-B",
                   lower = c(1e-5, 1e-5), upper = c(Inf, Inf))
      gamma_mle = fit$par[1]
      Lambda_mle = fit$par[2]
      if( any(is.na(fit$par)) || any(fit$par < 0) ){
        U_FD = 1
      }else{
        U_FD <- exp(compute_log_UBMarkov_FD( RmaxFD, gamma_mle, Lambda_mle, Kobs_i, sum(Nj_i), alpha, M_max ))
        # U_FD = 1
      }
      U_FD <- min(1,U_FD)
      if (!is.na(U_FD) && U_FD < eps) {
        stopped_FD <- TRUE
        Nstop_FD   <- ni
      }
      # stopped_FD = TRUE
    }
    
    ## ---- 5.2 PYP (on observed data) ----
    if (!stopped_PYP) {
      start_params <- c(alpha = 0.1, theta = 1)
      fit <- optim(par = start_params, fn = llik_pyp,
                   n = sum(Nj_i), Kn = Kobs_i, data_obs = Nj_i, # extra parameters
                   method = "L-BFGS-B",
                   lower = c(0, -1), upper = c(1-1e-10, Inf))
      alpha_mle = fit$par[1]
      theta_mle = fit$par[2]
      if( any(is.na(fit$par)) || any(fit$par < 0) ){
        U_PYP = 1
      }else{
        U_PYP <- exp(compute_log_UBMarkov( Rmax, alpha_mle, theta_mle, Kobs_i, sum(Nj_i), alpha ))
      }
      U_PYP <- min(1,U_PYP)
      if (!is.na(U_PYP) && U_PYP <= eps) {
        stopped_PYP <- TRUE
        Nstop_PYP   <- ni
      }
      # stopped_PYP = TRUE
    }
    
    ## ---- 5.3 Freq (on observed data) ----
    if (!stopped_Freq) {
      U_Freq <- ub_pain(n = sum(Nj_i), Rmax = Rmax, alfa = alpha)
      U_Freq <- min(1,U_Freq)
      if (!is.na(U_Freq) && U_Freq <= eps) {
        stopped_Freq <- TRUE
        Nstop_Freq   <- ni
      }
    }
    
    # Early exit if all four rules have stopped
    if (stopped_FD && stopped_PYP && stopped_Freq) break
  }
  
  ## ------------------------------------------------------------
  ## Post-loop: handle rules that *never* stopped by n_max
  ## ------------------------------------------------------------
  if (!stopped_FD) {
    Nstop_FD <- n_max
  }
  if (!stopped_PYP) {
    Nstop_PYP <- n_max
  }
  if (!stopped_Freq) {
    Nstop_Freq <- n_max
  }
  
  return(c(Nstop_FD,Nstop_PYP,Nstop_Freq))
}


SRabu_grid = function( eps_grid, data, nstart,
                       Nrep, num_cores, seed0,
                       alpha = 0.05, M_max = 200)
{
  Lgrid = length(eps_grid) # grid length
  res_list = vector("list",Lgrid)
  res_list = lapply(res_list, function(x) {
    y = matrix(nrow = Nrep, ncol = 3)
    colnames(y) = c("FD","PYP","Freq")
    y
  }  )
  
  
  ## Parallel run (no prints allowed)
  cluster <- makeCluster(num_cores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)
  clusterExport(cluster, list("SRabu_grid_single_run"),
                envir = environment())
  res_list = parLapply( cl = cluster, x = eps_grid,
                        fun = SRabu_grid_multiple_run,
                        data = data, nstart = nstart,
                        alpha = alpha, M_max = M_max,
                        seed0 = seed0, Nrep = Nrep)
  stopCluster(cluster)
  
  return(res_list)
}



SRinc_grid_multiple_run <- function(eps, data, nstart, seed0, Nrep, alpha)
{
  ## Functions
  suppressWarnings(suppressPackageStartupMessages(library(tibble)))
  suppressWarnings(suppressPackageStartupMessages(library(parallel)))
  suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
  suppressWarnings(suppressPackageStartupMessages(library(progress)))
  suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
  source("../../R/Rfunctions.R")
  Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
  source("../../../BinomialCIs/R/Rfunctions.R")
  Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")
  
  
  set.seed(seed0)
  seeds = sample(1:999999, size = Nrep)
  res = matrix(nrow = Nrep, ncol = 5)
  colnames(res) = c("3IBP","MixPois","MixBin","Freq.Bdd","Freq.Ubd")
  
  for(ii in 1:Nrep){
    seed = seeds[ii]
    res[ii,] = SRinc_grid_single_run(eps, data, nstart, seed, alpha)
  }
  return(res)
}

SRinc_grid_single_run <- function(eps, data, nstart, seed, alpha)
{
  set.seed(seed) # set seed
  RmaxFD = 50; Rmax = 100; beta = 1e-5;
  var_gamma = 10; var_nb    = 10
  
  n = nrow(data) # total num. obs.
  Kn = ncol(data) # total num. distinct
  ordered_idx = sample(1:n, size = n) # choose ordering of obs.
  
  if(nstart >= (n-1))
    stop("nstart must be smaller than n-1")
  
  # Stopping flags and outputs
  stopped_3IBP <- stopped_MixPois <- stopped_MixBin <- stopped_FreqBdd <- stopped_FreqUbd <- FALSE
  Nstop_3IBP <- Nstop_MixPois <- Nstop_MixBin <- Nstop_FreqBdd <- Nstop_FreqUbd <- NA_integer_
  
  ## ------------------------------------------------------------
  ## Run loop up to n_max = n
  ## ------------------------------------------------------------
  n_max = n
  ni = 2
  for(ni in nstart:(n_max-1)) {
    # Allow for a non-multiple n_max if needed
    remaining <- n_max - ni
    if (remaining <= 0L) break
    
    ## ---- Observed vector (true + error species) ----
    idx_species_i = ordered_idx[1:ni] # select obs. up to time ni
    data_i = data[idx_species_i,] 
    Nj_i = colSums(data_i) # compute frequencies
    Nj_i = Nj_i[Nj_i > 0]
    Kobs_i = length(which(Nj_i > 0))
    if( Kobs_i == 0L) next   # nothing observed yet
    
    ## ---- 5.1 3IBP (on observed data) ----
    if (!stopped_3IBP) {
      # Param. estimation (3 params PP)
      start_params <- c(alpha = 0.1, gamma= 1, u = 1)
      fit <- optim(par = start_params, fn = llik_PP3Parm, 
                   method = "L-BFGS-B",
                   n = ni, Kn = Kobs_i, data_obs = Nj_i,
                   lower = c(1e-16, 1e-16, 1e-16), 
                   upper = c(1-1e-10, Inf, Inf)) 
      alpha_mle = fit$par[1]
      gamma_mle = fit$par[2]
      c_mle     = fit$par[3] - alpha_mle
      
      # Upper bound (3 params PP)
      if( any(is.na(fit$par)) || any(fit$par < 0) ){
        ub3IBP = 1
      }else{
        ub3IBP = exp(compute_log_UBMarkov_BeBePois(Rmax, alpha_mle, c_mle, gamma_mle, ni, alpha ))
        # ub3IBP = 1
      }
      ub3IBP <- min(1,ub3IBP)
      if (!is.na(ub3IBP) && ub3IBP < eps) {
        stopped_3IBP <- TRUE
        Nstop_3IBP   <- ni
      }
      # stopped_3IBP = TRUE
    }
    
    ## ---- 5.2 MixPois (on observed data) ----
    if (!stopped_MixPois) {
      # Param. estimation (3 params PP)
      start_params <- c(alpha = 0.1, u = 1, mu_gamma = 1)
      fit <- optim(par = start_params, fn = llik_MixPois,
                   method = "L-BFGS-B",
                   n = ni, Kn = Kobs_i, data_obs = Nj_i, var_gamma = var_gamma,
                   lower = c(1e-16, 1e-16, 1e-16), upper = c(1-1e-10, Inf, Inf))
      alpha_mle = fit$par[1]
      c_mle = fit$par[2] - alpha_mle
      mugamma_mle = fit$par[3]
      
      gamma_hyperparams = gamma_shape_rate(mugamma_mle,var_gamma)
      u = gamma_hyperparams$shape
      v = gamma_hyperparams$rate
      
      
      # Upper bound (3 params PP)
      if( any(is.na( c(alpha_mle, c_mle,u,v) )) || any(c(alpha_mle, c_mle,u,v) < 0) ){
        ubMixPois = 1
      }else{
        ubMixPois = exp(compute_log_UBMarkov_BeBeMixPois( Rmax, alpha_mle, c_mle, ni, Kobs_i, u, v, alpha))
      }
      ubMixPois <- min(1,ubMixPois)
      if (!is.na(ubMixPois) && ubMixPois < eps) {
        stopped_MixPois <- TRUE
        Nstop_MixPois   <- ni
      }
      # stopped_MixPois = TRUE
    }
    
    ## ---- 5.3 MixBin (on observed data) ----
    if (!stopped_MixBin) {
      # Param. estimation (3 params PP)
      start_params <- c(a = 1, b = 1, mu_nb = 1)
      fit <- optim(par = start_params, fn = llik_MixBin,
                   method = "L-BFGS-B",
                   n = ni, Kn = Kobs_i, data_obs = Nj_i, var_nb = var_nb,
                   lower = c(1e-10, 1e-10, 1e-10), upper = c(Inf, Inf, var_nb-1e-10))
      a_mle = fit$par[1]
      b_mle = fit$par[2] 
      munb_mle = fit$par[3]
      
      nb_hyperparams = NegBin_params(munb_mle,var_nb)
      r_nb = nb_hyperparams$r
      p_nb = nb_hyperparams$p
      
      if( any(is.na(c(a_mle,b_mle,r_nb,p_nb))) || any(c(a_mle,b_mle,r_nb,p_nb) < 0) ){
        ubMixBin = 1
      }else{
        ubMixBin = exp(compute_log_UBMarkov_BeBeMixNBin( Rmax, a_mle, b_mle, ni, Kobs_i, r_nb, p_nb, alpha))
        # U_FD = 1
      }
      ubMixBin <- min(1,ubMixBin)
      if (!is.na(ubMixBin) && ubMixBin < eps) {
        stopped_MixBin <- TRUE
        Nstop_MixBin   <- ni
      }
      # stopped_MixBin = TRUE
    }
    
    ## ---- 5.4 Freq.Bdd (on observed data) ----
    if (!stopped_FreqBdd) {
      b_n <- log(ni)
      Mguess = 10 * Kobs_i
      Nj_guess = c(Nj_i, rep(0,Mguess - length(Nj_i) ))
      ubFreqBdd <- compute_UB_analytical(ni, Nj_guess, Mguess, b_n, alpha, FALSE)
      ubFreqBdd = min(1,ubFreqBdd); ubFreqBdd = max(0,ubFreqBdd)
      if (!is.na(ubFreqBdd) && ubFreqBdd <= eps) {
        stopped_FreqBdd <- TRUE
        Nstop_FreqBdd   <- ni
      }
    }
    
    ## ---- 5.5 Freq.Ubd (on observed data) ----
    if (!stopped_FreqUbd) {
      Shat  <- sum(Nj_i) / ni
      Sstar <- ( sqrt( -log(beta) / (2 * ni) ) +
                   sqrt( Shat + (-log(beta) / (2 * ni)) ) )^2
      r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) + log(ni) - log(log(ni))
      ubFreqUbd <- compute_UB_rnorm(ni, alpha, beta, r_n, Shat)
      ubFreqUbd = min(1,ubFreqUbd); ubFreqUbd = max(0,ubFreqUbd)
      
      # xx = data.frame("eps" = eps,"ni" = ni, "ubFreqUbd" = ubFreqUbd)
      # trimmed_eps = get_first3digits(eps,3)
      # save(xx, file = paste0("temp/","eps",trimmed_eps,"_ni",ni,".Rdat"))
      
      if (!is.na(ubFreqUbd) && ubFreqUbd <= eps) {
        stopped_FreqUbd <- TRUE
        Nstop_FreqUbd   <- ni
      }
    }
    
    # Early exit if all four rules have stopped
    if (stopped_3IBP && stopped_MixPois && stopped_MixBin && stopped_FreqBdd && stopped_FreqUbd) break
  }
  
  ## ------------------------------------------------------------
  ## Post-loop: handle rules that *never* stopped by n_max
  ## ------------------------------------------------------------
  if (!stopped_3IBP) {
    stopped_3IBP <- n_max
  }
  if (!stopped_MixPois) {
    stopped_MixPois <- n_max
  }
  if (!stopped_MixBin) {
    stopped_MixBin <- n_max
  }
  if (!stopped_FreqBdd) {
    stopped_FreqBdd <- n_max
  }
  if (!stopped_FreqUbd) {
    stopped_FreqUbd <- n_max
  }
  
  return( c(Nstop_3IBP,Nstop_MixPois,Nstop_MixBin,Nstop_FreqBdd,Nstop_FreqUbd) )
}

SRinc_grid = function( eps_grid, data, nstart,
                       Nrep, num_cores, seed0,
                       alpha = 0.05)
{
  Lgrid = length(eps_grid) # grid length
  res_list = vector("list",Lgrid)
  res_list = lapply(res_list, function(x) {
    y = matrix(nrow = Nrep, ncol = 5)
    colnames(y) = c("3IBP","MixPois","MixBin","Freq.Bdd","Freq.Ubd")
    y
  }  )
  
  
  ## Parallel run (no prints allowed)
  cluster <- makeCluster(num_cores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)
  clusterCall(cluster, getwd)
  clusterExport(cluster, list("SRinc_grid_single_run"),
                envir = environment())
  res_list = parLapply( cl = cluster, x = eps_grid,
                        fun = SRinc_grid_multiple_run,
                        data = data, nstart = nstart,
                        alpha = alpha,
                        seed0 = seed0, Nrep = Nrep)
  stopCluster(cluster)
  
  return(res_list)
}

# Coverage based stopping rules -------------------------------------------

## Coverages stopping rules
SRabu_cov_grid_multiple_run <- function(cov, data, nstart, seed0, Nrep)
{
  ## Functions
  suppressWarnings(suppressPackageStartupMessages(library(tibble)))
  suppressWarnings(suppressPackageStartupMessages(library(parallel)))
  suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
  suppressWarnings(suppressPackageStartupMessages(library(progress)))
  suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
  source("../../R/Rfunctions.R")
  Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
  
  set.seed(seed0)
  seeds = sample(1:999999, size = Nrep)
  res = matrix(nrow = Nrep, ncol = 1)
  colnames(res) = c("Coverage")
  
  for(ii in 1:Nrep){
    seed = seeds[ii]
    res[ii,] = SRabu_cov_grid_single_run(cov, data, nstart, seed)
  }
  return(res)
}

SRabu_cov_grid_single_run <- function(cov, data, nstart ,seed)
{
  set.seed(seed) # set seed
  
  nj = table(data) # compute abundances
  n  = sum(nj) # total num. obs.
  Kn = length(nj) # total num. distinct
  ordered_idx = sample(1:n, size = n) # choose ordering of obs.
  
  # Stopping flags and outputs
  stopped_cov <- FALSE
  Nstop_cov   <- NA_integer_
  
  ## ------------------------------------------------------------
  ## Run loop up to n_max = n
  ## ------------------------------------------------------------
  n_max = n
  ni = 2
  for(ni in nstart:(n_max-1)) {
    # Allow for a non-multiple n_max if needed
    remaining <- n_max - ni
    if (remaining <= 0L) break
    
    ## ---- Observed abundance vector (true + error species) ----
    idx_species_i = ordered_idx[1:ni] # select obs. up to time ni
    data_i = data[idx_species_i] 
    Nj_i = table(data_i) # compute frequencies
    Nj_i = Nj_i[Nj_i > 0]
    Kobs_i = length(which(Nj_i > 0))
    if( Kobs_i == 0L) next   # nothing observed yet
    
    ## ---- Coverage-based rule ----
    if (!stopped_cov) {
      C_hat <- SpadeR:::Chat.Ind(Nj_i, m = sum(Nj_i))   
      if (!is.na(C_hat) && C_hat >= cov) {
        stopped_cov <- TRUE
        Nstop_cov   <- ni
      }
    }
    
    # Early exit if all four rules have stopped
    if (stopped_cov) break
  }
  
  ## ------------------------------------------------------------
  ## Post-loop: handle rules that *never* stopped by n_max
  ## ------------------------------------------------------------
  if (!stopped_cov) {
    Nstop_cov <- n_max
  }
  
  return(c(Nstop_cov))
}


SRabu_cov_grid = function(cov_grid, data, nstart, Nrep, num_cores, seed0)
{
  Lgrid = length(cov_grid) # grid length
  res_list = vector("list",Lgrid)
  res_list = lapply(res_list, function(x) {
    y = matrix(nrow = Nrep, ncol = 1)
    colnames(y) = c("Coverage")
    y
  }  )
  
  
  ## Parallel run (no prints allowed)
  cluster <- makeCluster(num_cores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)
  clusterExport(cluster, list("SRabu_cov_grid_single_run"),
                envir = environment())
  res_list = parLapply( cl = cluster, x = cov_grid,
                        fun = SRabu_cov_grid_multiple_run,
                        data = data, nstart = nstart,
                        seed0 = seed0, Nrep = Nrep)
  stopCluster(cluster)
  
  return(res_list)
}

# Lorenzo Ghilotti's functions for extrapolation curves -------------------------------------
convert_features_list <- function(feature_matrix){
  
  feat_list <- vector("list", nrow(feature_matrix))
  
  for (i in 1:nrow(feature_matrix)){
    feat_list[[i]] <- which(feature_matrix[i,]==1, arr.ind = TRUE)
  }
  
  return (feat_list)
}
rarefaction.array <- function(object, n_reorderings = 1, seed = 1234) {
  
  feature_list <- convert_features_list(object)
  n <- nrow(object)
  
  if (n_reorderings == 1){
    
    rare_curve <- sapply(1:n, function(i) length(unique(unlist(feature_list[1:i]))) )
    
  } else {
    
    rare_curve <- matrix(NA, nrow = n_reorderings, ncol = n)
    
    for (j in 1:n_reorderings){
      
      f_list <- sample(feature_list)
      
      rare_curve[j, ] <- sapply(1:n, function(i) length(unique(unlist(f_list[1:i]))) )
    }
    
    # rare_curve <- colMeans(rare_curve)
    rare_curve_qnt = apply(rare_curve,2,quantile, prob = c(0.025,0.5,0.975))
  }
  
  return(rare_curve_qnt)
  
}



# Utilities ---------------------------------------------------------------
gamma_moments = function(a,b){
  list("mean" = a/b, "var" = a/(b*b))
}
gamma_shape_rate = function(mu,sig2){
  list("shape" = (mu*mu)/sig2, "rate" = mu/(sig2))
}

NegBin_moments = function(r,p){
  list("mean" = r*(1-p)/(p), "var" = r*(1-p)/(p*p))
}
NegBin_params = function(mu,sig2){
  list("p" = (mu)/sig2, "r" = (mu*mu)/(sig2-mu))
}
