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

library("ProductFormFA")

llik_BBNBin = function(x, n, Kn, data_obs, one_over_p) {
  log_a <- x[1]; log_b <- x[2]; # Beta hyperparameters
  logmu_nb <- x[3];             # Neg. Bin. hyperparameters
  a <- exp(log_a); b <- exp(log_b);
  mu_nb <- exp(logmu_nb); p_nb <- 1/one_over_p
  r_nb <- mu_nb * (1-p_nb)/p_nb  
  -log_efpfBeBeMixNBin( n, Kn, data_obs, a, b, r_nb, p_nb )
}
neg_llik_BBNBin = function(x, n, Kn, data_obs, one_over_p) {
  a <- x[1]; b <- x[2]; # Beta hyperparameters
  mu_nb <- x[3];             # Neg. Bin. hyperparameters
  # a <- exp(log_a); b <- exp(log_b);
  # mu_nb <- exp(logmu_nb); 
  p_nb <- 1/one_over_p
  r_nb <- mu_nb * (1-p_nb)/p_nb  
  if(mu_nb > 1e5)
    -exp(20)
  log_efpfBeBeMixNBin( n, Kn, data_obs, a, b, r_nb, p_nb )
}

my_GibbsFA_eb <- function(feature_matrix, model, type, seed = 1234, 
                          eb_params = NULL, 
                          Nhat_MM = NULL, var_fct = NULL, 
                          var_GammaIBP = NULL, ...) 
{
  
  if (type == "EFPF"){
    
    cat("\n Handle data ... ")
    # Remove NAs and 0s column from feature_matrix
    feature_matrix <- feature_matrix[, colSums(is.na(feature_matrix))==0]
    feature_matrix <- feature_matrix[, colSums(feature_matrix)!=0]
    counts <- colSums(feature_matrix)
    n <- nrow(feature_matrix)
    K <- ncol(feature_matrix)
    cat("done\n")
    if (all(class(eb_params) == c("eb_params", "BB"))){ # optimizing alpha, theta, Nhat_prime
      cat("\n Beta Bernoulli case \n")
      if (model == "PoissonBB_eb" | model == "NegBinBB_eb" | model == "classicBB_eb") {
        
        eb_init <- eb_params$init # this contains Nhat_prime, NOT Nhat
        eb_known <- eb_params$known
        cat("\n Optimization ... ")
        res <- nlminb(
          start = eb_init, objective =  my_neg_log_EFPF_BB_R,  
          n = n, counts = counts, known = eb_known, 
          lower = c(-Inf, 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf)
        )
        cat("done\n")
        Nhat_res <- round(unname(res$par["Nhat_prime"]) + K)
        cat("\n Nhat_res = ",Nhat_res,"\n")
        if (model == "PoissonBB_eb"){
          cat("\n Beta Bernoulli case --> PoissonBB_eb \n")
          out <- list("feature_matrix" = feature_matrix,
                      "eb_params" = eb_params,
                      "alpha" = unname(res$par["alpha"]), 
                      "theta" = unname(res$par["s"] - res$par["alpha"]),
                      "lambda" = Nhat_res,
                      "fun_value" = res$objective
          )
          
          class(out) <- c("GibbsFA", "PoissonBB_eb")
          return(out)
        }
        
        if (model == "NegBinBB_eb"){
          cat("\n Beta Bernoulli case --> NegBinBB_eb \n")
          out <- list("feature_matrix" = feature_matrix,
                      "eb_params" = eb_params,
                      "alpha" = unname(res$par["alpha"]), 
                      "theta" = unname(res$par["s"] - res$par["alpha"]),
                      "var_fct" = var_fct,
                      "n0" = Nhat_res/(var_fct - 1),
                      "mu0" = Nhat_res,
                      "fun_value" = res$objective
          )
          
          class(out) <- c("GibbsFA", "NegBinBB_eb")
          return(out)
        }
        
        if (model == "classicBB_eb"){
          cat("\n Beta Bernoulli case --> classicBB_eb \n")
          out <- list("feature_matrix" = feature_matrix,
                      "eb_params" = eb_params,
                      "alpha" = unname(res$par["alpha"]), 
                      "theta" = unname(res$par["s"] - res$par["alpha"]),
                      "N" = Nhat_res,
                      "fun_value" = res$objective
          )
          
          class(out) <- c("GibbsFA", "classicBB_eb")
          return(out)
        }
      }
      cat("\n Beta Bernoulli case ---> DONE! \n")
    } 
    
    if (all(class(eb_params) == c("eb_params", "IBP"))){ # optimizing alpha, theta, Gamma
      cat("\n IBP case \n")
      if (model == "GammaIBP_eb" | model == "classicIBP_eb"){ # we always optimize all the parameters, never fix them
        
        eb_init <- eb_params$init 
        eb_known <- eb_params$known
        cat("\n Optimization ... ")
        res <- nlminb(
          start = eb_init, objective =  neg_log_EFPF_IBP_R,  
          n = n, counts = counts, known = eb_known, 
          lower = c(1e-5, 1e-5, 1e-5), upper = c(1 - 1e-5, Inf, Inf)
        )
        cat("done\n")
        Gamma_prior_mean <- unname(res$par["Gamma"])
        cat("\n Gamma_prior_mean = ",Gamma_prior_mean,"\n")
        if (model == "GammaIBP_eb"){
          cat("\n IBP case --> GammaIBP_eb \n")
          out <- list("feature_matrix" = feature_matrix,
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
          cat("\n IBP case --> classicIBP_eb \n")
          out <- list("feature_matrix" = feature_matrix,
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
      cat("\n IBP case ---> DONE!! \n")
    }
    
    
    if (model == "PoissonBB") {
      cat("\n Model PoissonBB \n")
      # Initialization of the optimization
      eb_init <- eb_params$init
      eb_known <- eb_params$known
      cat("\n Optimization ... ")
      res <- nlminb(
        start = eb_init, objective =  neg_log_EFPF_GibbsFA_R, model = "PoissonBB",
        n = n, counts = counts, known = eb_known, lower = c(-Inf, 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf)
      )
      cat(" Done! \n")
      
      out <- list("feature_matrix" = feature_matrix,
                  "eb_params" = eb_params,
                  "alpha" = unname(res$par["alpha"]),
                  "theta" = unname(res$par["s"] - res$par["alpha"]),
                  "lambda" = eb_known[["lambda"]],
                  "fun_value" = res$objective
      )
      
      class(out) <- c("GibbsFA", "PoissonBB_eb")
      return(out)
    }
    
    if (model == "NegBinBB") {
      cat("\n Model NegBinBB \n")
      # Initialization of the optimization
      eb_init <- eb_params$init
      eb_known <- eb_params$known
      cat("\n Optimization ... ")
      res <- nlminb(
        start = eb_init, objective = neg_log_EFPF_GibbsFA_R, model = "NegBinBB",
        n = n, counts = counts, known = eb_known,
        lower = c(-Inf, 1e-5, 1 + 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf, Inf)
      )
      cat(" Done! \n")
      
      out <- list("feature_matrix" = feature_matrix,
                  "eb_params" = eb_params,
                  "alpha" = unname(res$par["alpha"]),
                  "theta" = unname(res$par["s"] - res$par["alpha"]),
                  "var_fct" = eb_known[["var_fct"]],
                  "n0" = eb_known[["mu0"]]/(eb_known[["var_fct"]] - 1),
                  "mu0" = eb_known[["mu0"]],
                  "fun_value" = res$objective
      )
      
      class(out) <- c("GibbsFA", "NegBinBB_eb")
      return(out)
    }
    
    if (model == "GammaIBP") {
      cat("\n Model NegBinBB \n")
      stop("not implemented")
    }
    
    
    stop("Incompatible eb_params and model parameters.")
    
  }
}

my_neg_log_EFPF_BB_R <- function(n, counts, par, known){ # par: alpha, s, Nhat' = Nhat - k
  
  par <- ifelse(is.na(known), par, known)
  
  return(neg_log_EFPF_BB(n, counts, par))
  
}
# Check -------------------------------------------------------------------

# Beta parameters
a = 0.5; b = 10
alpha = -a; s = b; theta = a+b

# Neg. Binomial parameters
mu = 120;v  = 240
r_nb = (mu*mu)/(v-mu)
p_nb = mu/v
q_nb = 1-p_nb

params_PF = c(alpha,s,r_nb,p_nb)

M_mc = rnbinom(n = 10000, size = r_nb, prob = p_nb)
MM = max(M_mc)
ptrue = rbeta(n=M_mc[1],shape1 = a, shape2 = b)

# Gen. data
n = 1000
Nj = rbinom(n = length(ptrue), size = n, prob = ptrue) 
Kn = length(which(Nj > 0))
N_j = Nj[which(Nj > 0)]


log_efpfBeBeMixNBin( n, Kn, N_j, a, b, r_nb, p_nb )
neg_log_EFPF_NegBinBB(n,N_j,params_PF)



# Optimization -------------------------------------------------------------------
seed = 34231
set.seed(seed)

## Gen. data -------------------------------------------------------------------
Mb = length(ptrue)
data <- matrix(rbinom(n * Mb, 1, ptrue), nrow = n, ncol = Mb, byrow = TRUE)
data = data[,colSums(data) > 0]
n = nrow(data)
Kn = ncol(data)
N_j = colSums(data)
names(N_j) = as.character(1:Kn)
Nj_ordered = sort(N_j, decreasing = TRUE)



## PFFA -------------------------------------------------------------------
eb_init_BB <- list(alpha = -1, s = 100, Nhat_prime = 50)
eb_known_BB <- list()
eb_params_obj_BB <- eb_params(model = "BB", init = eb_init_BB, known = eb_known_BB )

res = my_GibbsFA_eb(feature_matrix = data,
                    model = "NegBinBB_eb", type = "EFPF",
                    eb_params =  eb_params_obj_BB, 
                    var_fct = 10)

res$fun_value

aest = -res$alpha
best = res$theta - aest
n0est = res$n0
mu0est = res$mu0
pest = n0est/(n0est+mu0est)
1/pest

round(c(aest,best,mu0est),4)

M_mc = rnbinom(n = 10000, size = n0est, prob = pest)
xmax = 50
w_mc_list = lapply(M_mc, function(m) {w = rbeta(n = m, shape1 = aest, shape2 = best); w = sort(w,decreasing = TRUE); w[1:xmax]})
mat_MixBin = do.call(cbind, w_mc_list)
mat_MixBin[is.na(mat_MixBin)] = 0
qnt_MixBin = apply(mat_MixBin, 1, quantile, probs = c(0.025,0.5,0.975), na.rm = TRUE)

ymax = max(max(Nj_ordered)/n , max(qnt_MixBin) )
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_MixBin[1,1:xmax], rev(qnt_MixBin[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:Kn, y = Nj_ordered/n, col = "darkred", pch = 16, cex = 0.5)

## Ale -------------------------------------------------------------------

one_over_p = 10
start_params <- c( log_a = log(1), log_b = log(1), logmu_nb = log(Kn) )
fit <- optim(par = start_params, fn = llik_BBNBin,
             method = "L-BFGS-B",
             n = n, Kn = Kn, data_obs = N_j, one_over_p = one_over_p, 
             lower = c(-Inf,-Inf, -Inf), upper = c(Inf, Inf, Inf))
fit$value



a_mle = exp(fit$par[1])
b_mle = exp(fit$par[2])
mu_mle = exp(fit$par[3])
p_mle = 1/one_over_p
r_mle = mu_mle * (1-p_mle)/p_mle

round(c(a_mle,b_mle,mu_mle),4)


M_mc = rnbinom(n = 10000, size = r_mle, prob = p_mle)
xmax = 50
w_mc_list = lapply(M_mc, function(m) {w = rbeta(n = m, shape1 = a_mle, shape2 = b_mle); w = sort(w,decreasing = TRUE); w[1:xmax]})
mat_MixBin = do.call(cbind, w_mc_list)
mat_MixBin[is.na(mat_MixBin)] = 0
qnt_MixBin = apply(mat_MixBin, 1, quantile, probs = c(0.025,0.5,0.975), na.rm = TRUE)

ymax = max(max(Nj_ordered)/n , max(qnt_MixBin) )
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2,2.5,1,0), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(qnt_MixBin[1,1:xmax], rev(qnt_MixBin[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:Kn, y = Nj_ordered/n, col = "darkred", pch = 16, cex = 0.5)




# Criminals --------------------------------------------------------------------

load("RawDataInc.Rdat")
n = ncol(A)
Kn = nrow(A)
N_j = rowSums(A)
names(N_j) = as.character(1:Kn)
Nj_ordered = sort(N_j, decreasing = TRUE)

seed = 34231
set.seed(seed)

data = t(A)

var_fct_vec = c(1,5,10,50,100,500,1000)
eb_init_BB <- list(alpha = -1, s = 100, Nhat_prime = 50)
eb_known_BB <- list()
eb_params_obj_BB <- eb_params(model = "BB", init = eb_init_BB, known = eb_known_BB )
res_vec_PF = lapply(var_fct_vec, function(x){
  fit <- tryCatch(GibbsFA_eb(feature_matrix = data,
                             model = "NegBinBB_eb", type = "EFPF",
                             eb_params =  eb_params_obj_BB, 
                             var_fct = x),
                  error=function(e) NA)
  res = list("a_mle" = -fit$alpha,
             "b_mle" = fit$theta+fit$alpha,
             "mu_mle" = fit$mu0,
             "p_mle" = 1/x,
             "r_mle" = fit$n0,
             "fval" = fit$fun_value)
  res
})

sapply(res_vec_PF, function(x) x$fval)
sapply(res_vec_PF, function(x) x$r_mle)
sapply(res_vec_PF, function(x) x$mu_mle)
sapply(res_vec_PF, function(x) x$a_mle)
sapply(res_vec_PF, function(x) x$b_mle)



one_over_p_vec = c(1,5,10,50,100,500,1000)
Ltrials = length(one_over_p_vec)

res_vec = lapply(one_over_p_vec, function(x){
  start_params <- c( log_a = log(1), log_b = log(1), logmu_nb = log(Kn) )
  fit <- tryCatch(optim(par = start_params, fn = llik_BBNBin,
                        method = "L-BFGS-B",
                        n = n, Kn = Kn, data_obs = N_j, one_over_p = x, 
                        lower = c(-Inf,-Inf, -Inf), upper = c(Inf, Inf, Inf)),
                  error=function(e) NA)
  res = list("a_mle" = exp(fit$par[1]),
             "b_mle" = exp(fit$par[2]),
             "mu_mle" = exp(fit$par[3]),
             "p_mle" = 1/x,
             "r_mle" = exp(fit$par[3]) * (1-(1/x))/(1/x),
             "fval" = fit$value)
  res
})

sapply(res_vec, function(x) x$fval)
sapply(res_vec, function(x) x$r_mle)
sapply(res_vec, function(x) x$mu_mle)
sapply(res_vec, function(x) x$a_mle)
sapply(res_vec, function(x) x$b_mle)



