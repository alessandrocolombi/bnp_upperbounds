# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C)  Alessandro Colombi (alessandro.colombi@unibocconi.it)                        |
#                                                                                                 |
#     Part of this code is adapted from ProductFormFA package by Lorenzo Ghilotti                 |
#     downloaded from https://github.com/LGhilotti/ProductFormFA/tree/main)                       |
#                                                                                                 |
#     This repo is free software: you can redistribute it and/or modify it under                  |
#     the terms of the GNU General Public License as published by the Free                        |
#     Software Foundation; see <https:#cran.r-project.org/web/licenses/GPL-3>.                    |
#                                                                                                 |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
library("ProductFormFA")
my_neg_log_EFPF_BB_R <- function(n, counts, par, known){ # par: alpha, s, Nhat' = Nhat - k
  
  par <- ifelse(is.na(known), par, known)
  
  return(neg_log_EFPF_BB(n, counts, par))
  
}
my_neg_log_EFPF_IBP_R <- function(n, counts, par, known){ # par: alpha, s, Gamma
  
  par <- ifelse(is.na(known), par, known)
  
  return(neg_log_EFPF_IBP(n, counts, par))
  
}

opt_GibbsFA_eb <- function(n, counts, model, type, 
                           eb_params = NULL, 
                           Nhat_MM = NULL, var_fct = NULL, 
                           var_GammaIBP = NULL, ...) 
{
  
  if (type == "EFPF"){
    
    # Remove NAs and 0s from counts
    counts <- counts[which(!is.na(counts))]
    counts <- counts[which(counts > 0)]
    K <- length(counts)
    if (all(class(eb_params) == c("eb_params", "BB"))){ # optimizing alpha, theta, Nhat_prime
      
      if (model == "PoissonBB_eb" | model == "NegBinBB_eb" | model == "classicBB_eb") {
        
        eb_init <- eb_params$init # this contains Nhat_prime, NOT Nhat
        eb_known <- eb_params$known
        res <- nlminb( start = eb_init, objective =  my_neg_log_EFPF_BB_R,  
                       n = n, counts = counts, known = eb_known, 
                       lower = c(-Inf, 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf) )
        Nhat_res <- round(unname(res$par["Nhat_prime"]) + K)
        if (model == "PoissonBB_eb"){
          out <- list("counts" = counts,
                      "eb_params" = eb_params,
                      "alpha" = unname(res$par["alpha"]), 
                      "theta" = unname(res$par["s"] - res$par["alpha"]),
                      "lambda" = Nhat_res,
                      "fun_value" = res$objective )
          
          class(out) <- c("GibbsFA", "PoissonBB_eb")
          return(out)
        }
        
        if (model == "NegBinBB_eb"){
          out <- list("counts" = counts,
                      "eb_params" = eb_params,
                      "alpha" = unname(res$par["alpha"]), 
                      "theta" = unname(res$par["s"] - res$par["alpha"]),
                      "var_fct" = var_fct,
                      "n0" = Nhat_res/(var_fct - 1),
                      "mu0" = Nhat_res,
                      "fun_value" = res$objective )
          
          class(out) <- c("GibbsFA", "NegBinBB_eb")
          return(out)
        }
        
        if (model == "classicBB_eb"){
          out <- list("counts" = counts,
                      "eb_params" = eb_params,
                      "alpha" = unname(res$par["alpha"]), 
                      "theta" = unname(res$par["s"] - res$par["alpha"]),
                      "N" = Nhat_res,
                      "fun_value" = res$objective)
          
          class(out) <- c("GibbsFA", "classicBB_eb")
          return(out)
        }
      }
    } 
    
    if (all(class(eb_params) == c("eb_params", "IBP"))){ # optimizing alpha, theta, Gamma
      if (model == "GammaIBP_eb" | model == "classicIBP_eb"){ # we always optimize all the parameters, never fix them
        
        eb_init <- eb_params$init 
        eb_known <- eb_params$known
        res <- nlminb( start = eb_init, objective =  my_neg_log_EFPF_IBP_R,  
                       n = n, counts = counts, known = eb_known, 
                       lower = c(1e-5, 1e-5, 1e-5), upper = c(1 - 1e-5, Inf, Inf))
        Gamma_prior_mean <- unname(res$par["Gamma"])
        if (model == "GammaIBP_eb"){
          out <- list("counts" = counts,
                      "eb_params" = eb_params,
                      "alpha" = unname(res$par["alpha"]), 
                      "theta" = unname(res$par["s"] - res$par["alpha"]),
                      "var" = var_GammaIBP,
                      "gam" = Gamma_prior_mean,
                      "a" = Gamma_prior_mean^2 / var_GammaIBP,
                      "b" = Gamma_prior_mean / var_GammaIBP,
                      "fun_value" = res$objective
          )
          
          class(out) <- c("GibbsFA", "GammaIBP_eb")
          return(out)
        }
        
        if (model == "classicIBP_eb"){
          out <- list("counts" = counts,
                      "eb_params" = eb_params,
                      "alpha" = unname(res$par["alpha"]), 
                      "theta" = unname(res$par["s"] - res$par["alpha"]),
                      "gam" = Gamma_prior_mean,
                      "fun_value" = res$objective
          )
          
          class(out) <- c("GibbsFA", "classicIBP_eb")
          return(out)
        }
      }
    }
    
    stop("Incompatible eb_params and model parameters.")
    
  }else{
    stop("type must be equal to EFPF")
  }
}


ParEst_PFFA <- function(n,counts,model,
                        alpha0=NULL,s0=NULL,Nhat0=NULL,Gamma0=NULL,
                        Nhat_MM = NULL, var_fct = NULL, 
                        var_GammaIBP = NULL, ...)
{
  eb_known_BB <- eb_known_IBP <- list()
  if (model == "PoissonBB_eb" | model == "NegBinBB_eb" | model == "classicBB_eb"){
    # Set initial parameters
    if(is.null(alpha0))
      alpha0 = -1
    if(is.null(s0))
      s0 = 100
    if(is.null(Nhat0))
      Nhat0 = 50
    eb_init_BB <- list(alpha = alpha0, s = s0, Nhat_prime = Nhat0)
    
    # cat("\n ",eb_init_BB$alpha," ",eb_init_BB$s," ",eb_init_BB$Nhat_prime)
    # Define seach object
    eb_params_obj_BB <- eb_params(model = "BB", init = eb_init_BB, known = eb_known_BB )
    
    out = opt_GibbsFA_eb(n = n, counts = counts,
                         model = model, type = "EFPF",
                         eb_params =  eb_params_obj_BB, 
                         var_fct = var_fct)
    
    a_mle = -out$alpha
    b_mle = out$theta - a_mle
    r_nb_mle = out$n0
    mu_nb_mle = out$mu0
    p_nb_mle = r_nb_mle/(r_nb_mle+mu_nb_mle)
    fval = out$fun_value
    res = list("a_mle" = a_mle,
               "b_mle" = b_mle,
               "r_nb_mle" = r_nb_mle,
               "mu_nb_mle" = mu_nb_mle,
               "p_nb_mle" = p_nb_mle,
               "fval" = fval)
    return(res)
  }
  if (model == "GammaIBP_eb" | model == "classicIBP_eb"){
    # Set initial parameters
    if(is.null(alpha0))
      alpha0 = 0.5
    if(is.null(s0))
      s0 = 1
    if(is.null(Gamma0))
      Gamma0 = 10
    eb_init_IBP <- list(alpha = alpha0, s = s0, Gamma = Gamma0)
    
    # Define seach object
    eb_params_obj_IBP <- eb_params(model = "IBP", init = eb_init_IBP, known = eb_known_IBP )
    
    out = opt_GibbsFA_eb(n = n, counts = counts,
                         model = model, type = "EFPF",
                         eb_params =  eb_params_obj_IBP, 
                         var_GammaIBP = var_GammaIBP)
    sigma_mle = out$alpha
    c_mle = out$theta 
    gam_mle = out$gam
    fval = out$fun_value
    res = list("sigma_mle" = sigma_mle, "c_mle" = c_mle,
               "gam_mle" = gam_mle, "fval" = fval)
    return(res)
  }
}















