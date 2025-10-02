# wd = "C:/Users/colom/bnp_upperbounds/Rscripts/CancerData/"
# setwd(wd)
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

# Load --------------------------------------------------------------------
load("cancer_types.Rdat")
load("cancer_names_easy.Rdat")

# idx = 2

cancer_name = cancer_names[idx]
filename = paste0(wd,"TCGA/",cancer_types[idx],"_targeted.RData")
load(filename)
cat("\n idx = ",idx,"; cancer name: ",cancer_name,"\n")

seed = 42
set.seed(seed)

exp_name = paste0("Cancer_",cancer_types[idx])
save_exp = TRUE
save_img_1 = TRUE
save_img_2 = TRUE
file_name = paste0("save/",exp_name,".Rdat")
img_name_1 = paste0("img/",exp_name,"_len.pdf")
img_name_2 = paste0("img/",exp_name,"_nlen.pdf")

# Read & Statistics -------------------------------------------------------

Nj = c(apply(Z, 2, sum))
n = nrow(Z)
Kobs = length( which(Nj > 0) )
TabNj = c(length(which(Nj == 0)), tabulate(Nj, nbins = n) )



# Options -----------------------------------------------------------------
train_prop = seq(0.05,1,by = 0.05)
Nexp = length(train_prop)

Nrep_max = 200
Nrep = ceiling(seq(Nrep_max,1,length.out = Nexp))

alfa = 0.05
beta = 1e-5
Rmax = 100

# Run ---------------------------------------------------------------------

# Save structure for 3IBP
ub_3IBP_list <- lapply(1:Nexp, function(ii){rep(NA,Nrep[ii])})
nub_3IBP_list <- ub_3IBP_list

# Save structure for Freq. estimators
ub_Inter_list <- lapply(1:Nexp, function(ii){rep(NA,Nrep[ii])})
ub_rnorm_list <- nub_Inter_list <- nub_rnorm_list <- ub_Inter_list


for(ii in 1:Nexp){ # Loop for each training size
  cat("\n ii = ",ii,"\n")
  prop = train_prop[ii]
  n_train = ceiling( n*prop )
  for(b in 1:Nrep[ii]){ # For each train, repeat Nrep[ii] times
    
    idx_train = sample( 1:n, size = n_train ) # select rows to include in train
    idx_train = sort(idx_train) 
    Z_train = Z[idx_train,] # find training set
    Nj_train = c(apply(Z_train, 2, sum)) # summarise training set
    Kobs_train = length( which(Nj_train > 0) ) # num. distinct features in training
    data_obs = Nj_train[which(Nj_train > 0)] # species frequencies
    

    
    #1.1) Param. estimation (3 params PP)
    start_params <- c(alpha = 0.1, gamma= 1, u = 1)
    fit <- optim(par = start_params, fn = llik_PP3Parm, 
                 method = "L-BFGS-B",
                 n = n_train, Kn = Kobs_train, data_obs = data_obs,
                 lower = c(1e-16, 1e-16, 1e-16), 
                 upper = c(1-1e-10, Inf, Inf)) 
    alpha_mle = fit$par[1]
    gamma_mle = fit$par[2]
    c_mle     = fit$par[3] - alpha_mle
    
    ##1.2) Upper bound (3 params PP)
    ub_3IBP_list[[ii]][b] = exp(compute_log_UBMarkov_BeBePois( Rmax, alpha_mle, c_mle, gamma_mle, n_train, alfa ))
    nub_3IBP_list[[ii]][b] = n_train * ub_3IBP_list[[ii]][b]
    
    #2.1) Compute Shat
    Shat = sum( Nj_train )/n_train
    
    ##2.2) Intersection
    ub_Inter_list[[ii]][b] = compute_UB_intersection(n_train, alfa, beta, Shat)
    nub_Inter_list[[ii]][b] = n_train * ub_Inter_list[[ii]][b]
    
    ##2.3) r-norm
    Sstar = ( sqrt( -log(beta)/(2*n_train) ) + sqrt( Shat + (-log(beta)/(2*n_train)) ) )^2
    rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n_train) - log(log(n_train))
    ub_rnorm_list[[ii]][b] = compute_UB_rnorm(n_train, alfa, beta, rn, Shat)
    nub_rnorm_list[[ii]][b] = n_train * ub_rnorm_list[[ii]][b]
  }
}

res = list( "ub_3IBP_list" = ub_3IBP_list, "nub_3IBP_list" = nub_3IBP_list,
            "ub_Inter_list" = ub_Inter_list, "nub_Inter_list" = nub_Inter_list,
            "ub_rnorm_list" = ub_rnorm_list, "nub_rnorm_list" = nub_rnorm_list,
            "n_train_vec" = ceiling( n*train_prop ) )


if(save_exp)
  save(res, file = file_name)



# Final summary and plot ---------------------------------------------------

# 3IBP)
ub_3IBP_list  = lapply(ub_3IBP_list, quantile, probs = c(0.025,0.5,0.975)); ub_3IBP_list = do.call(cbind, ub_3IBP_list)
nub_3IBP_list  = lapply(nub_3IBP_list, quantile, probs = c(0.025,0.5,0.975)); nub_3IBP_list = do.call(cbind, nub_3IBP_list)

# Freq. est.) 
ub_Inter  = lapply(ub_Inter_list, quantile, probs = c(0.025,0.5,0.975)); ub_Inter = do.call(cbind, ub_Inter)
ub_rnorm  = lapply(ub_rnorm_list, quantile, probs = c(0.025,0.5,0.975)); ub_rnorm = do.call(cbind, ub_rnorm)
nub_Inter = lapply(nub_Inter_list, quantile, probs = c(0.025,0.5,0.975)); nub_Inter = do.call(cbind, nub_Inter)
nub_rnorm = lapply(nub_rnorm_list, quantile, probs = c(0.025,0.5,0.975)); nub_rnorm = do.call(cbind, nub_rnorm)

## plot (1) : CI length ----------------------------------------------------

ymax = max( c( ub_3IBP_list[3,],ub_Inter[3,],ub_rnorm[3,]) ) * 1.05
ymin = 0 #min( c( nub_Inter[1,],nub_rnorm[1,])) 
ylabs = round(seq(ymin,ymax,length.out = 5),1)

n_train_vec = ceiling( n*train_prop )



if(save_img_1)
  pdf(img_name_1)
par( mfrow = c(1,1), mar = c(3.5,3.5,1,0.5), mgp=c(2.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "n", ylab = "Length",
     xlim = c(0,max(n_train_vec) ) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = n_train_vec, y = ub_3IBP_list[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkblue" )
polygon( c(n_train_vec, rev(n_train_vec)),
         c(ub_3IBP_list[1,], rev(ub_3IBP_list[3,])),
         col = ACutils::t_col("darkblue",30),
         border = NA) # plot in-sample bands
points( x = n_train_vec, y = ub_Inter[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkred" )
polygon( c(n_train_vec, rev(n_train_vec)),
         c(ub_Inter[1,], rev(ub_Inter[3,])),
         col = ACutils::t_col("darkred",30),
         border = NA) # plot in-sample bands
points( x = n_train_vec, y = ub_rnorm[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkgreen" ) 
polygon( c(n_train_vec, rev(n_train_vec)),
         c(ub_rnorm[1,], rev(ub_rnorm[3,])),
         col = ACutils::t_col("darkgreen",30),
         border = NA) # plot in-sample bands
legend("topright",c("3IBP","Intersection","r-norm"), 
       lwd = 3, col = c("darkblue","darkred","darkgreen"))
if(save_img_1)
  dev.off()


## plot (2) : CI length x n ----------------------------------------------------

ymax = max( c( nub_3IBP_list[3,] ) ) * 1.05
ymin = 0 #min( c( nub_Inter[1,],nub_rnorm[1,])) 
ylabs = round(seq(ymin,ymax,length.out = 5),1)

n_train_vec = ceiling( n*train_prop )


if(save_img_2)
  pdf(img_name_2)
par( mfrow = c(1,1), mar = c(3.5,3.5,1,0.5), mgp=c(2.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "n", ylab = "n x Length",
     xlim = c(0,max(n_train_vec) ) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = n_train_vec, y = nub_3IBP_list[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkblue" )
polygon( c(n_train_vec, rev(n_train_vec)),
         c(nub_3IBP_list[1,], rev(nub_3IBP_list[3,])),
         col = ACutils::t_col("darkblue",30),
         border = NA) # plot in-sample bands
points( x = n_train_vec, y = nub_Inter[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkred" )
polygon( c(n_train_vec, rev(n_train_vec)),
         c(nub_Inter[1,], rev(nub_Inter[3,])),
         col = ACutils::t_col("darkred",30),
         border = NA) # plot in-sample bands
points( x = n_train_vec, y = nub_rnorm[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkgreen" ) 
polygon( c(n_train_vec, rev(n_train_vec)),
         c(nub_rnorm[1,], rev(nub_rnorm[3,])),
         col = ACutils::t_col("darkgreen",30),
         border = NA) # plot in-sample bands
legend("topright",c("3IBP","Intersection","r-norm"), 
       lwd = 3, col = c("darkblue","darkred","darkgreen"))
if(save_img_2)
  dev.off()





