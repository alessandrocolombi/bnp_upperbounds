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

eps_test = 0.15
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
