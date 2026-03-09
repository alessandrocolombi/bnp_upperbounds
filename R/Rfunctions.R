
# Gen. Dist. (species) ----------------------------------------------------
sim_zipfs = function(M,s){
  w = sapply(1:M,function(j) j^(-s))
  w / sum(w)
}
sim_geom = function(M,a){
  w = sapply(1:M,function(j) (a)^(j)  )
  w / sum(w)
}
sim_negbin = function(M,l,r){
  w = dnbinom(x = 0:(M-1), size = l, prob = r) # Painsky parametrization
  w / sum(w)
}
sim_betabin = function(M,a,b){
  w = dbetabinom.ab(x = 0:(M-1), size = M-1, shape1 = a, shape2 = b)
  w / sum(w)
}
sim_unif = function(M){
  rep(1/M,M)
}


# Worst case distribution -------------------------------------------------

find_ma_worstunif <- function(n,alpha){
  lStir = lastirlings2(n)
  for(m in 1:n){
    lres = lfactorial(m) - n*log(m) + lStir[(n+1),(m+1)]
    prob = 1 - exp(lres)
    if(prob > alpha){
      break
    }
  }
  m
}
worst_uniform <- function(M,n,alpha){
  ma = find_ma_worstunif(n,alpha)
  if(M <= ma)
    p = rep(1/M,M)
  if(M > ma){
    p = rep(0,M)
    p[1:ma] = rep(1/ma, ma)
  }
  p
}

p_all_seen_uniform <- function(n, ma) { 
  if (ma < 0) return(0) 
  if (n < ma) return(0) # impossible to see all symbols 
  lStir = lastirlings2(n) 
  lres = lfactorial(ma) - n*log(ma) + lStir[(n+1),(ma+1)] 
  exp(lres) 
} 
oracle_worst_uniform <- function(n, ma, alpha) {
  1/ma
  # q <- 1 - p_all_seen_uniform(n, ma) # P(missing at least one) 
  # if(q > alpha) 1/ma else 0 
}


# General sim from species  -----------------------------------------------

sim_generic_species = function(name,M,params){
  if(name == "Zipfs"){
    sim_zipfs(M,params[1])
  }else if( name == "Geom" ){
    sim_geom(M,params[1])
  }else if( name == "Uniform"){
    sim_unif(M)
  }
  else if( name == "NegBin"){
    sim_negbin(M, params[1], params[2])
  }
  else if( name == "BetaBin"){
    sim_betabin(M,params[1], params[2])
  }
  else if( name == "WorstUnif"){
    worst_uniform(M,params[1],params[2])
  }
  else
    stop("Invalid name")
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
llik_DirMult = function(x, n, M, data){
  gamma = x[1]
  -log_DirMulti(n,M,data,gamma)
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
    Nstop_3IBP <- n_max
  }
  if (!stopped_MixPois) {
    Nstop_MixPois <- n_max
  }
  if (!stopped_MixBin) {
    Nstop_MixBin <- n_max
  }
  if (!stopped_FreqBdd) {
    Nstop_FreqBdd <- n_max
  }
  if (!stopped_FreqUbd) {
    Nstop_FreqUbd <- n_max
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


# Ferguson-Klass for beta process -----------------------------------------

base_rnd <- function(n) runif(n)

#' Ferguson-Klass sampler for the 3-parameter (stable) Beta process
#'
#' Parameters:
#'   c     > 0            concentration / strength
#'   sigma in [0,1)       discount
#'   gamma > 0            
#'   base_rnd(n)          function sampling n i.i.d. atom locations from normalized B0
#'   eps                  stop when weights fall below eps (approximation)
#'   n_grid               grid size for precomputing tail integral
#'
#' Returns:
#'   list(w = weights in decreasing order, theta = sampled locations)
fk_stable_beta_process <- function(c, sigma, gamma, base_rnd = base_rnd,
                                   eps = 1e-8,
                                   n_grid = 20000L,
                                   grid_min = 1e-12,
                                   grid_max = 1 - 1e-12,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # ---- checks ----
  stopifnot(is.numeric(c), length(c) == 1L, is.finite(c), c > 0)
  stopifnot(is.numeric(sigma), length(sigma) == 1L, is.finite(sigma), sigma >= 0, sigma < 1)
  stopifnot(is.numeric(gamma), length(gamma) == 1L, is.finite(gamma), gamma > 0)
  stopifnot(is.numeric(eps), length(eps) == 1L, is.finite(eps), eps > 0, eps < 1)
  stopifnot(is.function(base_rnd))
  
  # ---- precompute tail integral T(x) on a grid ----
  # Use a grid dense near 0 because of the w^{-1-sigma} singularity.
  # log-spaced in [grid_min, 1), then clamp to grid_max.
  grid <- exp(seq(log(grid_min), log(grid_max), length.out = n_grid))
  grid[grid >= grid_max] <- grid_max
  
  # log Lévy density (up to constant):
  # f(w) = w^{-1-sigma} (1-w)^{c+sigma-1}
  logf <- (-(1 + sigma)) * log(grid) + (c + sigma - 1) * log1p(-grid)
  f <- exp(logf)
  
  # Tail integral via trapezoid rule from w to 1:
  # tail_raw[i] ~ ???_{grid[i]}^1 f(w) dw
  # compute cumulative integral from the right
  dx <- diff(grid)
  trap <- 0.5 * (f[-length(f)] + f[-1L]) * dx
  
  tail_raw <- numeric(length(grid))
  tail_raw[length(grid)] <- 0
  # tail at grid[i] includes trap[i] + ... + trap[end]
  tail_raw[-length(grid)] <- rev(cumsum(rev(trap)))
  
  # Multiply by constant gamma  
  const <- gamma #* c / gamma(1 - sigma)
  tail <- const * tail_raw
  
  # We need inverse of tail: given t in (0, tail(grid_min)] find w with T(w)=t.
  # tail decreases with w; for approxfun x must be increasing.
  tail_inc <- rev(tail)   # increasing
  grid_dec <- rev(grid)   # corresponding w values increasing in tail_inc's index
  
  inv_tail <- approxfun(x = tail_inc, y = grid_dec, method = "linear",
                        rule = 2, ties = "ordered")
  
  # Stopping threshold in Gamma-space: stop when Gamma_i > T(eps)
  # If eps is smaller than grid_min, clamp to grid_min.
  eps_clamped <- max(eps, grid_min)
  # approximate T(eps) by interpolating tail on grid
  T_eps <- approx(x = grid, y = tail, xout = eps_clamped, rule = 2)$y
  
  if (!is.finite(T_eps) || T_eps <= 0) {
    return(list(w = numeric(0), theta = base_rnd(0)))
  }
  
  # ---- Ferguson-Klass series ----
  w_list <- numeric(0)
  G <- 0.0
  i <- 0L
  
  repeat {
    i <- i + 1L
    G <- G + rexp(1L, rate = 1)  # Gamma_i
    if (G > T_eps) break
    
    w_i <- inv_tail(G)
    # numerical safety
    if (!is.finite(w_i) || w_i <= 0 || w_i >= 1) next
    if (w_i < eps) break
    w_list <- c(w_list, w_i)
  }
  
  # Sample locations
  theta <- base_rnd(length(w_list))
  
  list(w = w_list, theta = theta)
}




# GOF ---------------------------------------------------------------------

SimModel_PD = function(Mmax,Nrep,params,seed = 091079){
  sigma = params[1]; theta = params[2]
  sim_PYP = r_SB(Nrep,Mmax,sigma,theta,seed)
  sim_PYP = apply(sim_PYP, 1, sort, decreasing = TRUE)
  # pyp_qnt = apply(sim_PYP, 1, quantile, probs = c(0.025,0.5,0.975))
  sim_PYP
}
SimModel_FDP = function(Mmax,Nrep,params,seed = 091079){
  set.seed(seed)
  gamma = params[1]; Lambda = params[2]
  M_mc = rpois(n=Nrep,lambda = Lambda) + 1
  sim_FD = matrix(0,nrow = Nrep, ncol = Mmax)
  ii = 1
  for(ii in 1:Nrep){
    w = rgamma(n = M_mc[ii], shape = gamma, rate = 1); w = w/sum(w); w = sort(w, decreasing = TRUE)
    len = min(Mmax,M_mc[ii])
    sim_FD[ii,1:len] = w[1:len]
  }
  sim_FD = apply(sim_FD, 1, sort, decreasing = TRUE)
  # FD_qnt = apply(sim_FD, 1, quantile, probs = c(0.025,0.5,0.975))
  sim_FD
}
SimModel_DirMulti = function(Mmax = NULL,Nrep,params,seed = 091079){
  set.seed(seed)
  gamma = params[1]; M = params[2]
  sim_DM = matrix(0,nrow = Nrep, ncol = M)
  ii = 1
  for(ii in 1:Nrep){
    w = rgamma(n = M, shape = gamma, rate = 1); w = w/sum(w); w = sort(w, decreasing = TRUE)
    sim_DM[ii,] = w
  }
  sim_DM = apply(sim_DM, 1, sort, decreasing = TRUE)
  # DM_qnt = apply(sim_DM, 1, quantile, probs = c(0.025,0.5,0.975))
  sim_DM
}
SimModel_generic = function(model,Mmax,Nrep,params,seed = 091079){
  if(model == "PD"){
    SimModel_PD(Mmax,Nrep,params,seed)
  }else if(model == "FDP"){
    SimModel_FDP(Mmax,Nrep,params,seed)
  }else if(model == "DirMulti"){
    SimModel_DirMulti(Mmax,Nrep,params,seed)
  }else{
    stop("model must be PD or FDP or DirMulti")
  }
}

SimData = function(n,ptrue_mat, seed = 091079 ){
  set.seed(seed)
  M = nrow(ptrue_mat); Nrep = ncol(ptrue_mat)
  data_mat = apply(ptrue_mat, 2, function(ptrue) sample(1:M, size = n, replace = TRUE, prob = ptrue)) # n x Nrep
  data_mat
}

GOF_generic = function(model,n,Mmax,Nrep,params,
                       ptrue_mat = NULL,seed = 091079, AccCrv_length = 20){
  res_names = c("Envelop_qnt","Freq.Rare","AccCrv")
  res = vector("list", length = length(res_names))
  names(res) = res_names
  
  if(model == "True"){
    if(is.null(ptrue_mat))
      stop("ptrue_mat must be a Nrep x M matrix if model is set to true")
  }else{
    # Simulate probs from the model
    ptrue_mat = SimModel_generic(model,Mmax,Nrep,params,seed)
  }
  # Generate Nrep datasets
  data_mat = SimData(n,ptrue_mat, seed)
  
  # a) Envelop plot
  res$Envelop_qnt = apply(ptrue_mat, 1, quantile, probs = c(0.025,0.5,0.975))
  # b) Frequency of rare species
  res$Freq.Rare = apply(data_mat, 2, function(data){
    n_i = tabulate(data, nbins = Mmax)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    data_obs = n_i[idx_obs]
    r = 5; fr = rep(0,r); names(fr) = as.character(1:r)
    for(i in 1:r){
      fr[i] = length(which(data_obs == i))
    }
    fr
  })
  # c) Accumulation curve
  if( n < 10)
    stop("n must be larger than 10")
  ngrid = round(seq(10,n,length.out = AccCrv_length))
  mat <- lapply(ngrid, function(nn) {
    data_mat_n <- SimData(nn,ptrue_mat, seed)
    Kn_all = apply(data_mat_n, 2, function(data_n){
      n_i = tabulate(data_n, nbins = Mmax)
      idx_obs = which(n_i > 0)
      length(idx_obs)
    })
  })
  mat <- do.call(rbind,mat) 
  AccCrv = apply(mat, 1, quantile, probs = c(0.025,0.5,0.975), na.rm = TRUE)
  res$AccCrv = AccCrv
  
  res
}



SimModel_Post_PD = function(Kn,Nj,params,Mmax,Nrep,seed = 091079){
  if(length(Nj) != Kn)
    stop("Nj must be a vector of length Kn, absolute frequencies are required")
  
  sigma_PD = params[1]; theta_PD = params[2]
  shapes_PD = c(Nj - sigma_PD, theta_PD + Kn*sigma_PD)
  betas_PD = matrix(0,nrow = Nrep, ncol = Kn+1)
  betas_PD = t(apply(betas_PD,1,function(x){
    y = rgamma(n = Kn+1, shape = shapes_PD, rate = 1); y/sum(y)
  }))
  Pprime = t(SimModel_PD(Mmax,Nrep,c(sigma_PD,theta_PD)))
  Prob_post_PD = t(apply(cbind(betas_PD,Pprime), 1, function(x){
    res = rep(0,length(x)-1)
    res[1:Kn] = x[1:Kn]
    beta_prime = x[Kn+1]
    temp = x[(Kn+2):length(x)]
    res[(Kn+1):length(res)] = beta_prime * temp
    res
  }))
  Prob_post_PD
}
SimModel_Post_FDP = function(Kn,Nj,params,Mmax,Nrep,seed = 091079){
  if(length(Nj) != Kn)
    stop("Nj must be a vector of length Kn, absolute frequencies are required")
  
  n = sum(Nj)
  gamma_FDP = params[1]; Lambda_FDP = params[2]
  log_PMstar = rep(-Inf,Mmax) # Whole distribution (log-scale)
  log_PMstar = sapply(0:(Mmax-1), function(m) log_qMpost(m,n,Kn,gamma_FDP,Lambda_FDP,500))
  Mstar_MC = sample(0:(Mmax-1), size = Nrep, prob = exp(log_PMstar), replace = TRUE)
  Prob_post_FDP = matrix(0,nrow = Nrep, ncol = Kn+Mmax)
  for(b in 1:Nrep){
    Mstar = Mstar_MC[b]
    wj_unnorm = c(Nj+gamma_FDP,rep(gamma_FDP,Mstar))
    w = rgamma(n = (Kn+Mstar), shape = wj_unnorm, rate = 1); w = w/sum(w)
    Prob_post_FDP[b,1:(Kn+Mstar)] = w
  }
  Prob_post_FDP
}
SimModel_Post_DM = function(Kn,Nj,params,Mmax = NULL,Nrep,seed = 091079){
  gamma_DM = params[1]; M = params[2]
  wj_unnorm <- Nj + gamma_DM
  gamma_mat <- matrix(
    rgamma(n = M * Nrep, shape = rep(wj_unnorm, times = Nrep), rate = 1),
    nrow = M, ncol = Nrep
  )
  gamma_mat <- sweep(gamma_mat, 2, colSums(gamma_mat), FUN = "/")
  Prob_post_DM <- t(apply(gamma_mat, 2, sort, decreasing = TRUE))
  # Prob_post_DM <- gamma_mat
  Prob_post_DM
}
SimModel_Post_generic = function(model,Kn,Nj,params,Mmax,Nrep,seed = 091079){
  if(model == "PD"){
    SimModel_Post_PD(Kn,Nj,params,Mmax,Nrep,seed )
  }else if(model == "FDP"){
    SimModel_Post_FDP(Kn,Nj,params,Mmax,Nrep,seed )
  }else if(model == "DirMulti"){
    SimModel_Post_DM(Kn,Nj,params,Mmax,Nrep,seed )
  }else{
    stop("model must be PD or FDP or DirMulti")
  }
}


# Expected Num Species ----------------------------------------------------

Expected_Kn = function(model, n, Kn, params, M=NULL, 
                       seed = 31231, M_max = 500, Kn_max = 400)
{
  if(model == "PD"){
    sigma = params[1]; theta = params[2]
    if(sigma < 1e-10){
      # Dirichlet process case
      ExpKn = theta * ( digamma(theta+n) - digamma(theta) )
    }else{
      # Pitman-Yor case
      logA = lgamma(theta) - lgamma(theta+sigma) + lgamma(theta+n+sigma) - lgamma(theta+n) 
      ExpKn = theta/sigma * expm1(logA)
    }
    return( ExpKn )
  }else if(model == "FDP"){
    gamma = params[1]; Lambda = params[2]
    logC_vec = compute_logC(n, -gamma, 0.0)
    logC_vec = logC_vec[-1] # length n
    logV_vec = compute_logV_all( n, gamma, Lambda, M_max, Kn_max ) # length Kn_max (<= n)
    probs = exp( logC_vec[1:length(logV_vec)] + logV_vec )
    ExpKn = sum( probs * 1:length(probs) )
    return(ExpKn)
  }else if(model == "DirMulti"){
    stop("Not yet implemented")
  }
  stop("model must either be PD, FDP or DirMulti")
}

# ParEst - MCMC -----------------------------------------------------------


ParEst_MCMC_generic = function(model, n, Kn, Nj, Niter,
                               init_val, hy_prior, Adp_var,
                               UpdateParam, M=NULL, 
                               seed = 31231, M_max = 500)
{
  if(model == "FDP"){
    # Initial values
    gamma0 = init_val[1];  Lambda0 = init_val[2]
    # Hyperprior params
    a_gamma = hy_prior[1]; b_gamma = hy_prior[2];
    a_Lambda = hy_prior[3];b_Lambda = hy_prior[4];
    # Adaptive variance
    AdpVar_gamma = Adp_var[1]; AdpVar_Lambda = Adp_var[2]
    # Update or not
    UpdateGamma=UpdateParam[1]; UpdateLambda = UpdateParam[2]
    # Run
    res = GibbsSampler_FDP(n,Kn,Nj,Niter,gamma0,Lambda0,a_gamma,b_gamma,a_Lambda,b_Lambda,
                           AdpVar_gamma,AdpVar_Lambda,UpdateGamma,UpdateLambda,M_max,seed)
    return(res)
  }
  if(model == "DirMulti"){
    if(is.null(M))
      stop("M must be provided if model is DirMulti")
    
    # Initial values
    gamma0 = init_val[1]
    # Hyperprior params
    a_gamma = hy_prior[1];b_gamma = hy_prior[2]
    # Adaptive variance
    AdpVar_gamma = Adp_var[1]
    # Update or not
    UpdateGamma=UpdateParam[1]
    # Run
    res = GibbsSampler_DirMulti(n,Nj,M,Niter,gamma0,a_gamma,b_gamma,
                                AdpVar_gamma,UpdateGamma,seed)
    return(res)
  }
  if(model == "PD"){
    
    # Initial values
    sigma0 = init_val[1];theta0 = init_val[2]
    # Hyperprior params
    a_sigma = hy_prior[1];b_sigma = hy_prior[2];
    a_theta = hy_prior[3];b_theta = hy_prior[4];
    # Adaptive variance
    AdpVar_sigma = Adp_var[1]; AdpVar_theta = Adp_var[2]
    # Update or not
    UpdateSigma=UpdateParam[1]; UpdateTheta = UpdateParam[2]
    # Run
    res = GibbsSampler_PYP(n,Kn,Nj,Niter,sigma0,theta0,a_sigma,b_sigma,a_theta,b_theta,
                           AdpVar_sigma,AdpVar_theta,UpdateSigma,UpdateTheta,seed)
    return(res)
  }
  stop("model must either be PD, FDP or DirMulti")
}



Hychoice_MCMC_general = function(model, n, Kn, 
                                 hy_prior, M=NULL, 
                                 seed = 31231, M_max = 500)
{
  if(model == "PD"){
    a_sigma = hy_prior[1];b_sigma = hy_prior[2];
    ExpSigma = a_sigma/(a_sigma+b_sigma)
    fmin = function(x){ (Kn - Expected_Kn(model,n,Kn,params=c(ExpSigma,x)))^2 }
    start_params <- c(x = 1)
    fit <- optim(par = start_params, fn = fmin, 
                 method = "L-BFGS-B",
                 lower = c(1e-10), upper = c(1e10)) 
    ExpTheta = fit$par
    return( ExpTheta )
  }else if(model == "FDP"){
    ExpLambda = hy_prior[1]
    fmin = function(x){ (Kn - Expected_Kn(model,n,Kn,params=c(x,ExpLambda)))^2 }
    start_params <- c(x = 1)
    fit <- optim(par = start_params, fn = fmin, 
                 method = "L-BFGS-B",
                 lower = c(1e-10), upper = c(10000)) 
    ExpGamma = fit$par
    return( ExpGamma )
  }else if(model == "DirMulti"){
    stop("Not yet implemented")
  }
  stop("model must either be PD, FDP or DirMulti")
}










