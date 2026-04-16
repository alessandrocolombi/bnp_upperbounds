# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/CriminalData")
setwd(wd)

# Functions ---------------------------------------------------------------
source("../../R/Rfunctions.R")
source("../../R/PFFAfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")


# Options -----------------------------------------------------------------
seed0 = 34231
alpha <- alfa <- 0.05
Rmax = 100
subsample_grid = seq(10,45,by=5) # <--- modify here
Nrep = 50 # <--- modify here
save_res = FALSE # <--- modify here


# Load data ---------------------------------------------------------------
load("RawDataInc.Rdat")
data_full = t(A) # rows = people, cols = features
n_full = nrow(data_full)
Kn_full = ncol(data_full)

if (any(subsample_grid > n_full)) {
  stop("All subsample sizes must be <= the full sample size.")
}

set.seed(seed0)
subsample_seeds = sample(1:999999, size = length(subsample_grid) * Nrep)


# FB upper bound under subsampling ----------------------------------------
counter = 1
FB_subsampling_res = lapply(subsample_grid, function(n_sub) {
  cat("\n n_sub = ",n_sub," \n")
  rep_list = lapply(seq_len(Nrep), function(rep_idx) {
    cat(" ",rep_idx," - ")
    set.seed(subsample_seeds[counter])
    idx_sub = sort(sample.int(n_full, size = n_sub, replace = FALSE))
    data_sub = data_full[idx_sub, , drop = FALSE]

    N_j = colSums(data_sub)
    idx_obs = which(N_j > 0)
    Kn = length(idx_obs)
    Nj_obs = N_j[idx_obs]

    Mguess = 200 #Kn + M_offset
    Nj_guess = c(N_j, rep(0, Mguess - length(N_j)))

    start_params <- c(a = 1, b = 1)
    fit <- optim(
      par = start_params,
      fn = llik_FB,
      method = "L-BFGS-B",
      n = n_sub,
      Kn = Kn,
      data_obs = Nj_guess,
      M = Mguess,
      lower = c(1e-10, 1e-10),
      upper = c(Inf, Inf)
    )

    a_FB = fit$par[1]
    b_FB = fit$par[2]
    ubFB = exp(compute_log_UBMarkov_FB(Rmax, a_FB, b_FB, n_sub, Kn, Mguess, alpha))
    ubFB = min(ubFB, 1)

    counter <<- counter + 1

    data.frame(
      n_sub = n_sub,
      rep = rep_idx,
      Kn = Kn,
      Mguess = Mguess,
      a_FB = a_FB,
      b_FB = b_FB,
      omega = a_FB/(a_FB+b_FB),
      eta = (a_FB+b_FB+1),
      ubFB = ubFB,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rep_list)
})
FB_subsampling_res = do.call(rbind, FB_subsampling_res)

summarize_FB_subsampling = function(x) {
  mean_cols = setdiff(names(x), c("n_sub", "rep"))
  mean_vals = lapply(x[mean_cols], mean)
  names(mean_vals) = paste0("mean_", mean_cols)

  data.frame(
    n_sub = x$n_sub[1],
    as.data.frame(mean_vals, check.names = FALSE),
    q025_ubFB = quantile(x$ubFB, 0.025),
    q500_ubFB = quantile(x$ubFB, 0.5),
    q975_ubFB = quantile(x$ubFB, 0.975),
    row.names = NULL,
    check.names = FALSE
  )
}

FB_subsampling_summary = do.call(
  rbind,
  lapply(split(FB_subsampling_res, FB_subsampling_res$n_sub), summarize_FB_subsampling)
)
rownames(FB_subsampling_summary) = NULL


cat("\n=== FB subsampling results ===\n")
print(FB_subsampling_summary)

if (save_res) {
  save(FB_subsampling_res, FB_subsampling_summary,
       file = "save/Mod3_Test_FB_subsampling.Rdat")
}


# SR - Test ---------------------------------------------------------------

SRinc_test_FB <- function(eps, data, nstart, var_fct, seed, alpha)
{
  set.seed(seed) # set seed
  Rmax = 100
  
  n = nrow(data) # total num. obs.
  Kn = ncol(data) # total num. distinct
  ordered_idx = sample(1:n, size = n) # choose ordering of obs.
  
  if(nstart >= (n-1))
    stop("nstart must be smaller than n-1")
  
  # Stopping flag and output
  stopped_FB <- FALSE
  Nstop_FB <- NA_integer_
  
  ###
  ## Run loop up to n_max = n
  ###
  n_max = n
  ni = 10
  for(ni in nstart:(n_max-1)) {
    cat("\n ni = ",ni," - ")
    # Allow for a non-multiple n_max if needed
    remaining <- n_max - ni
    if (remaining <= 0L) break
    
    ### Observed vector (true + error species) 
    idx_species_i = ordered_idx[1:ni] # select obs. up to time ni
    data_i = data[idx_species_i,] 
    Nj_i = colSums(data_i) # compute frequencies
    Nj_i = Nj_i[Nj_i > 0]
    Kobs_i = length(which(Nj_i > 0))
    if( Kobs_i == 0L) next   # nothing observed yet
    
    ### Finite Beta (on observed data) 
    if (!stopped_FB) {
      # Param. estimation (Finite Beta)
      Mguess = 200
      Nj_guess = c(Nj_i, rep(0,Mguess - length(Nj_i) ))
      start_params <- c(a = 1, b = 1)
      fit <- optim(par = start_params, fn = llik_FB,
                   method = "L-BFGS-B",
                   n = ni, Kn = Kobs_i, data_obs = Nj_guess, M=Mguess,
                   lower = c(1e-10, 1e-10), upper = c(Inf, Inf))
      a_FB = fit$par[1]; 
      b_FB = fit$par[2];
      cat("(a,b) = ",a_FB,",",b_FB," - ")
      # Upper bound (FB)
      if( any(is.na( c(a_FB, b_FB) )) || any(c(a_FB, b_FB) < 0) ){
        ub_FB = 1
      }else{
        ub_FB = exp(compute_log_UBMarkov_FB( Rmax, a_FB, b_FB, ni, Kobs_i, Mguess, alpha))
      }
      ub_FB <- min(1,ub_FB)
      cat("ub_FB = ",ub_FB,", stop? ",ub_FB < eps," \n ")
      if (!is.na(ub_FB) && ub_FB < eps) {
        stopped_FB <- TRUE
        Nstop_FB   <- ni
      }
    }

    # Early exit once FB has stopped
    if (stopped_FB) break
  }
  
  ###
  ## Post-loop: handle the case that never stopped by n_max
  ###
  if (!stopped_FB) {
    Nstop_FB <- n_max
  }
  
  return(Nstop_FB)
}


# Execute -----------------------------------------------------------------
set.seed(seed0)

eps_test = 0.2
nstart_test = 10
var_fct_test = 100
seed_test = 12345

Nstop_FB_test = SRinc_test_FB(
  eps = eps_test,
  data = data_full,
  nstart = nstart_test,
  var_fct = var_fct_test,
  seed = seed_test,
  alpha = alpha
)

cat("\n=== FB stopping-rule test ===\n")
cat("eps =", eps_test, "| nstart =", nstart_test, "| Nstop_FB =", Nstop_FB_test, "\n")


# SR - Test - Old functions ---------------------------------------------------------------


SRinc_grid_multiple_run <- function(eps, data, nstart, var_fct, seed0, Nrep, alpha)
{
  ## Functions
  suppressWarnings(suppressPackageStartupMessages(library(tibble)))
  suppressWarnings(suppressPackageStartupMessages(library(parallel)))
  suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
  suppressWarnings(suppressPackageStartupMessages(library(progress)))
  suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
  
  source("../../R/Rfunctions.R")
  source("../../R/PFFAfunctions.R")
  Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
  
  # From BinomialCIs
  source("../../../BinomialCIs/R/Rfunctions.R")
  Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")
  
  
  set.seed(seed0)
  seeds = sample(1:999999, size = Nrep)
  res = matrix(nrow = Nrep, ncol = 5)
  colnames(res) = c("IBP","MBP","FB","Freq.Bdd","Freq.Ubd")
  
  for(ii in 1:Nrep){
    seed = seeds[ii]
    res[ii,] = SRinc_grid_single_run(eps=eps, data=data, nstart=nstart, var_fct=var_fct, seed=seed, alpha=alpha)
  }
  return(res)
}

SRinc_grid_single_run <- function(eps, data, nstart, var_fct, seed, alpha)
{
  set.seed(seed) # set seed
  RmaxFD = 50; Rmax = 100; beta = 1e-5;
  
  n = nrow(data) # total num. obs.
  Kn = ncol(data) # total num. distinct
  ordered_idx = sample(1:n, size = n) # choose ordering of obs.
  
  if(nstart >= (n-1))
    stop("nstart must be smaller than n-1")
  
  # Stopping flags and outputs
  stopped_3IBP <- stopped_FB <- stopped_MixBin <- stopped_FreqBdd <- stopped_FreqUbd <- FALSE
  Nstop_3IBP <- Nstop_FB <- Nstop_MixBin <- Nstop_FreqBdd <- Nstop_FreqUbd <- NA_integer_
  
  ###
  ## Run loop up to n_max = n
  ###
  n_max = n
  ni = 10
  for(ni in nstart:(n_max-1)) {
    # Allow for a non-multiple n_max if needed
    remaining <- n_max - ni
    if (remaining <= 0L) break
    
    ### Observed vector (true + error species) 
    idx_species_i = ordered_idx[1:ni] # select obs. up to time ni
    data_i = data[idx_species_i,] 
    Nj_i = colSums(data_i) # compute frequencies
    Nj_i = Nj_i[Nj_i > 0]
    Kobs_i = length(which(Nj_i > 0))
    if( Kobs_i == 0L) next   # nothing observed yet
    
    ### 5.1 3IBP (on observed data) 
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
    
    ### 5.2 Finite Beta (on observed data) 
    if (!stopped_FB) {
      cat("\n FB ... ")
      # Param. estimation (Finite Beta)
      Mguess = 200
      Nj_guess = c(Nj_i, rep(0,Mguess - length(Nj_i) ))
      start_params <- c(a = 1, b = 1)
      fit <- optim(par = start_params, fn = llik_FB,
                   method = "L-BFGS-B",
                   n = ni, Kn = Kobs_i, data_obs = Nj_guess, M=Mguess,
                   lower = c(1e-10, 1e-10), upper = c(Inf, Inf))
      a_FB = fit$par[1]; 
      b_FB = fit$par[2]; 
      # Upper bound (FB)
      if( any(is.na( c(a_FB, b_FB) )) || any(c(a_FB, b_FB) < 0) ){
        ub_FB = 1
      }else{
        ub_FB = exp(compute_log_UBMarkov_FB( Rmax, a_FB, b_FB, ni, Kobs_i, Mguess, alpha))
        cat(" ub_FB = ",ub_FB," ... ")
      }
      ub_FB <- min(1,ub_FB)
      if (!is.na(ub_FB) && ub_FB < eps) {
        cat("stopped ... ", ni, "\n")
        stopped_FB <- TRUE
        Nstop_FB   <- ni
      }
    }
    
    ### 5.3 MixBin (on observed data) 
    if (!stopped_MixBin) {
      # cat("\n MBP ... ")
      # Param. estimation (Mixed Binomial)
      eb_init_BB <- list(alpha = -1, s = 100, Nhat_prime = 500-Kobs_i)
      eb_known_BB <- list()
      eb_params_obj_BB <- eb_params(model = "BB", init = eb_init_BB, known = eb_known_BB )
      
      res = GibbsFA_eb(feature_matrix = data_i,
                       model = "NegBinBB_eb", type = "EFPF",
                       eb_params =  eb_params_obj_BB, 
                       var_fct = var_fct)
      
      a_mle = res$alpha+1
      b_mle = res$theta - a_mle
      r_nb = res$n0
      q_nb = 1 - 1/res$var_fct; p_nb = 1/var_fct
      
      if( any(is.na(c(a_mle,b_mle,r_nb,p_nb))) || any(c(a_mle,b_mle,r_nb,p_nb) < 0) ){
        ubMixBin = 1
      }else{
        ubMixBin = exp(compute_log_UBMarkov_BeBeMixNBin( Rmax, a_mle, b_mle, ni, Kobs_i, r_nb, p_nb, alpha))
        # cat(" ubMixBin = ",ubMixBin," ... ")
      }
      ubMixBin <- min(1,ubMixBin)
      if (!is.na(ubMixBin) && ubMixBin < eps) {
        # cat("stopped ... ", ni, "\n")
        stopped_MixBin <- TRUE
        Nstop_MixBin   <- ni
      }
      # stopped_MixBin = TRUE
    }
    
    ### 5.4 Freq.Bdd (on observed data) 
    if (!stopped_FreqBdd) {
      b_n <- log(ni)
      Mguess = 170 #10 * Kobs_i
      Nj_guess = c(Nj_i, rep(0,Mguess - length(Nj_i) ))
      ubFreqBdd <- compute_UB_analytical(ni, Nj_guess, Mguess, b_n, alpha, FALSE)
      ubFreqBdd = min(1,ubFreqBdd); ubFreqBdd = max(0,ubFreqBdd)
      if (!is.na(ubFreqBdd) && ubFreqBdd <= eps) {
        stopped_FreqBdd <- TRUE
        Nstop_FreqBdd   <- ni
      }
    }
    
    ### 5.5 Freq.Ubd (on observed data) 
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
    if (stopped_3IBP && stopped_FB && stopped_MixBin && stopped_FreqBdd && stopped_FreqUbd) break
  }
  
  ###
  ## Post-loop: handle rules that *never* stopped by n_max
  ###
  if (!stopped_3IBP) {
    Nstop_3IBP <- n_max
  }
  if (!stopped_FB) {
    Nstop_FB <- n_max
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
  
  return( c(Nstop_3IBP,Nstop_MixBin,Nstop_FB,
            Nstop_FreqBdd,Nstop_FreqUbd) )
}

SRinc_grid = function( eps_grid, data, nstart,
                       Nrep, num_cores, seed0,
                       var_fct = 100,
                       alpha = 0.05)
{
  Lgrid = length(eps_grid) # grid length
  res_list = vector("list",Lgrid)
  res_list = lapply(res_list, function(x) {
    y = matrix(nrow = Nrep, ncol = 5)
    colnames(y) = c("IBP","MBP","FB","Freq.Bdd","Freq.Ubd")
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
                        alpha = alpha, var_fct=var_fct,
                        seed0 = seed0, Nrep = Nrep)
  stopCluster(cluster)
  
  return(res_list)
}





# Load --------------------------------------------------------------------

load("RawDataInc.Rdat")
data = t(A)
n = ncol(A)
Kn = nrow(A)
N_j = rowSums(A)
names(N_j) = as.character(1:Kn)
Nj_ordered = sort(N_j, decreasing = TRUE)

seed = 34231
set.seed(seed)


# Options  --------------------------------------------------------
eps_grid = c(0.001, seq(0.1,0.3,length.out =  (34*5-1)) )
alpha = 0.05
M_max = 200
nstart = 10

seed0 = 4224
num_cores = 33 # <--- modify here
Nrep = 10 # <--- modify here

var_fct = 100
# Run) Mmax-based  --------------------------------------------------------
cat("\n Running stopping rule ... ")
res = SRinc_grid(eps_grid=eps_grid, data=data, nstart=nstart,
                 Nrep=Nrep, num_cores=num_cores, seed0=seed0,
                 alpha=alpha, var_fct=var_fct)
cat("done! Save and conclude \n")

