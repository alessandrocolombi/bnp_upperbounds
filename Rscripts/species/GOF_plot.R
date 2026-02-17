# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/species")
setwd(wd)

# Functions ---------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))

source("../../R/Rfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")


# Custom functions --------------------------------------------------------

# Plot options ------------------------------------------------------------------
save_img = FALSE
width = 12; height = 6
cex.labels <- cex.lab <- 2
cex.axis <- 2
cex.legend <- 1.5
mycol = c("darkorange","darkred","darkblue","lightpink","aquamarine")
mycol2 = c("black","lightblue")
lgd_names = c("Freq","PD","FDP","DirMulti-m","DirMulti")

# GOF ---------------------------------------------------------------------
params_zipfs = list(0.9,1.02,2)
params_geom = list(0.85,0.9,0.95)
params_unif = list(NA)
params_negbin = list(c(1,0.003),c(5,0.003),c(1,0.01))
experiments = list("Zipfs" = params_zipfs,
                   "Geom" = params_geom,
                   "Uniform" = params_unif,
                   "NegBin" = params_negbin)
Nrep = 200 # <---
n = 500
M = 200 
M_max = 200
seed = 232131
set.seed(seed)
Mmax = M; Nrep = 50
AccCrv_length = 50
ii = 4; name = names(experiments)[[ii]]
jj = 1
params = experiments[[ii]][[jj]]
trim_params = sapply(params, get_first3digits, 4)
if(length(trim_params)>1)
  trim_params = paste0(trim_params[1],"_",trim_params[2])


ptrue = sim_generic_species(name,M,params)
ptrue = sort(ptrue, decreasing = TRUE)

Mplot = 100
par( mfrow = c(1,1), mar = c(3.5,4.25,1,1), mgp=c(2.75,1,0), bty = "l", las = 1, cex.lab = 1 )
plot(x = 1:Mplot, y = ptrue[1:Mplot], type = "p", pch = 16, cex = 0.5,
     main = paste0(name,", ",ii,", ",jj), xlab = "", ylab = "")
ptrue_mat = matrix(nrow = Nrep,ncol = M)
ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)

# generate a single dataset
data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
n_i = tabulate(data, nbins = M)
idx_obs = which(n_i > 0)
Kn = length(idx_obs)
data_obs = n_i[idx_obs]

#0) TRUE
model = "True"
gof_true = GOF_generic(model,n,Mmax,Nrep,params,ptrue_mat = ptrue_mat, seed, AccCrv_length)

#a) PD
# Param. estimation (PYP)
start_params <- c(alpha = 0.5, theta = 1)
fit <- optim(par = start_params, fn = llik_pyp, 
             n = n, Kn = Kn, data_obs = data_obs, # extra parameters
             method = "L-BFGS-B",
             lower = c(0, -1), upper = c(1-1e-10, Inf)) 
alpha_mle = fit$par[1]
theta_mle = fit$par[2]
model = "PD"
params = c(alpha_mle,theta_mle)
gof_PD = GOF_generic(model=model,n=n,Mmax=Mmax,Nrep=Nrep,params=params,seed=seed,AccCrv_length=AccCrv_length)

#b) FDP
# Param. estimation (FD)
start_params <- c(gamma = 0.1, Lambda = Kn)
fit <- optim(par = start_params, fn = llik_FD, 
             n = n, Kn = Kn, data_obs = data_obs, M_max = M_max,# extra parameters
             method = "L-BFGS-B",
             lower = c(1e-5, 1e-5), upper = c(Inf, Inf)) 
gamma_mle = fit$par[1]
Lambda_mle = fit$par[2]
model = "FDP"
params = c(gamma_mle,Lambda_mle)
gof_FDP = GOF_generic(model=model,n=n,Mmax=Mmax,Nrep=Nrep,params=params,seed=seed,AccCrv_length=AccCrv_length)

#c) Dirichlet-Multinomial --> M = M
# Param. estimation (DirMulti)
start_params <- c(gamma = 0.5)
fit <- optim(par = start_params, fn = llik_DirMult,
             n = n, M = M, data = n_i, # extra parameters
             method = "L-BFGS-B",
             lower = c(1e-10), upper = c(Inf))
gamma_mle = fit$par[1]
params = c(gamma_mle,M)
gof_DM = GOF_generic(model=model,n=n,Mmax=Mmax,Nrep=Nrep,params=params,seed=seed,AccCrv_length=AccCrv_length)

# Plots

## a - Envelop
xmax = M
# ymax = max( c(ptrue, max(gof_PD$Envelop_qnt), max(gof_FDP$Envelop_qnt), max(gof_DM$Envelop_qnt)) )
ymax = max( ptrue ) * 1.33
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos*1e3, 0)

img_name = paste0("img/GOF/",name,"_",trim_params,"_GOFEnv",".pdf")
if(save_img)
  pdf(img_name, width = width, height = height)
par( mfrow = c(1,3), mar = c(3.5,6,1,1), mgp=c(4.4,1,0), bty = "l", las = 1, cex.lab = cex.lab )

# PD
plot(0,0,main = " ", ylab = "Prob * 1000 (PD)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n",
     xlab = "",cex.axis = cex.axis)
axis( side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis) 
polygon( c(1:xmax, rev(1:xmax)),
         c(gof_PD$Envelop_qnt[1,1:xmax], rev(gof_PD$Envelop_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:xmax, y = ptrue[1:xmax], col = "darkred", pch = 16)

# FDP
plot(0,0,main = " ", ylab = "Prob * 1000 (FDP)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n",xlab = "",cex.axis = cex.axis)
axis( side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis) 
polygon( c(1:xmax, rev(1:xmax)),
         c(gof_FDP$Envelop_qnt[1,1:xmax], rev(gof_FDP$Envelop_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:xmax, y = ptrue[1:xmax], col = "darkred", pch = 16)

# DM
plot(0,0,main = " ", ylab = "Prob * 1000 (DiriMulti)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n",xlab = "",cex.axis = cex.axis)
axis( side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis) 
polygon( c(1:xmax, rev(1:xmax)),
         c(gof_DM$Envelop_qnt[1,1:xmax], rev(gof_DM$Envelop_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:xmax, y = ptrue[1:xmax], col = "darkred", pch = 16)


if(save_img)
  dev.off()


## b - Frequency rare species
pos = c(4,8,12,16,20)
shift = c(-0.5,0,0.5)
xpos = pos; xlabs = as.character(1:5)
ymax = max(c( max(gof_true$Freq.Rare), max(gof_PD$Freq.Rare), max(gof_FDP$Freq.Rare) ))
ymin = 0
ylim_plot = c(ymin,ymax)

barplot_pos = c()
for(l in 1:length(pos)){
  barplot_pos = c(barplot_pos, pos[l]+shift)
}

img_name = paste0("img/GOF/",name,"_",trim_params,"_GOFRare",".pdf")
if(save_img)
  pdf(img_name, width = width, height = height)
par( mfrow = c(1,1), mar = c(3.5,4.5,1,1), mgp=c(3,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(0,22), ylim = ylim_plot,
     xlab = "r", ylab = paste0("Freq.Rare"), 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
grid(lty = 1,lwd = 1, col = "gray90" )
for(i in 1:nrow(gof_true$Freq.Rare)){
  boxplot(gof_true$Freq.Rare[i,], at = pos[i]+shift[2], add = T, 
          col = "black", pch = 16, yaxt = "n", cex = 0.5)
}
for(i in 1:nrow(gof_PD$Freq.Rare)){
  boxplot(gof_PD$Freq.Rare[i,], at = pos[i]+shift[1], add = T, 
          col = "darkred", pch = 16, yaxt = "n", cex = 0.5)
}
for(i in 1:nrow(gof_FDP$Freq.Rare)){
  boxplot(gof_FDP$Freq.Rare[i,], at = pos[i]+shift[3], add = T, 
          col = "darkblue", pch = 16, yaxt = "n", cex = 0.5)
}
legend("topright",c("True","PD","FDP"),
       fill = c("black","darkred","darkblue"), 
       cex = cex.legend, bty = "n", border = NA)
if(save_img)
  dev.off()


## c - Accumulation curve
ngrid = round(seq(10,n,length.out = AccCrv_length))
xpos = ngrid; xlabs = as.character(ngrid)
ymax = max(c( max(gof_true$AccCrv), max(gof_PD$AccCrv), max(gof_FDP$AccCrv) ))
ymin = 0
ylim_plot = c(ymin,ymax)


img_name = paste0("img/GOF/",name,"_",trim_params,"_GOFACC",".pdf")
if(save_img)
  pdf(img_name, width = width, height = height)
par( mfrow = c(1,1), mar = c(4,4.5,1,1), mgp=c(3,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(0,max(ngrid)), ylim = ylim_plot,
     xlab = "n", ylab = paste0(""), 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
polygon( c(ngrid, rev(ngrid)),
         c(gof_true$AccCrv[1,], rev(gof_true$AccCrv[3,])),
         col = scales::alpha("grey15", 0.25),
         border = NA) # plot in-sample bands
polygon( c(ngrid, rev(ngrid)),
         c(gof_PD$AccCrv[1,], rev(gof_PD$AccCrv[3,])),
         col = scales::alpha("darkred", 0.25),
         border = NA) # plot in-sample bands
polygon( c(ngrid, rev(ngrid)),
         c(gof_FDP$AccCrv[1,], rev(gof_FDP$AccCrv[3,])),
         col = scales::alpha("darkblue", 0.25),
         border = NA) # plot in-sample bands
points(x = ngrid, y = gof_true$AccCrv[2,], type = "l", lwd = 5, col = "black")
points(x = ngrid, y = gof_PD$AccCrv[2,], type = "l", lwd = 5, col = "darkred")
points(x = ngrid, y = gof_FDP$AccCrv[2,], type = "l", lwd = 5, col = "darkblue")
if(save_img)
  dev.off()




# Deviance ----------------------------------------------------------------

#a) Generate multiple datasets
data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
Deviance = t(apply(data_mat, 2, function(data){
  n_i = tabulate(data, nbins = M)
  idx_obs = which(n_i > 0)
  Kn = length(idx_obs)
  data_obs = n_i[idx_obs]
  
  res = matrix(nrow = 1, ncol = 3)
  colnames(res) = c("PD","FDP","DirMulti")
  #a) PD
  # Param. estimation (PYP)
  start_params <- c(alpha = 0.5, theta = 1)
  fit <- optim(par = start_params, fn = llik_pyp, 
               n = n, Kn = Kn, data_obs = data_obs, # extra parameters
               method = "L-BFGS-B",
               lower = c(0, -1), upper = c(1-1e-10, Inf)) 
  res[1,1] = 2*fit$value
  
  #b) FDP
  # Param. estimation (FDP)
  start_params <- c(gamma = 0.1, Lambda = Kn)
  fit <- optim(par = start_params, fn = llik_FD, 
               n = n, Kn = Kn, data_obs = data_obs, M_max = M_max,# extra parameters
               method = "L-BFGS-B",
               lower = c(1e-5, 1e-5), upper = c(Inf, Inf)) 
  gamma_mle = fit$par[1]
  Lambda_mle = fit$par[2]
  res[1,2] = 2*fit$value  
  
  #c) Dirichlet-Multinomial --> M = M
  # Param. estimation (DirMulti)
  start_params <- c(gamma = 0.5)
  fit <- optim(par = start_params, fn = llik_DirMult,
               n = n, M = M, data = n_i, # extra parameters
               method = "L-BFGS-B",
               lower = c(1e-10), upper = c(Inf))
  res[1,3] = 2*fit$value  
  
  res
})) # Nrep x 3 matrix

## b - Deviance plot
pos = c(2,4,6)
shift = c(0)
xpos = pos; xlabs = as.character(1:3)
ymax = max(Deviance);ymin = min(Deviance)
ylim_plot = c(ymin,ymax)

barplot_pos = c()
for(l in 1:length(pos)){
  barplot_pos = c(barplot_pos, pos[l]+shift)
}

img_name = paste0("img/GOF/",name,"_",trim_params,"_GOFDev",".pdf")
if(save_img)
  pdf(img_name, width = width, height = height)
par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(1,7), ylim = ylim_plot,
     xlab = "r", ylab = paste0("Deviance"), 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
grid(lty = 1,lwd = 1, col = "gray90" )
boxplot(Deviance[,1], at = pos[1], add = T, 
        col = "darkred", pch = 16, yaxt = "n", cex = 0.5)
boxplot(Deviance[,2], at = pos[2], add = T, 
        col = "darkblue", pch = 16, yaxt = "n", cex = 0.5)
boxplot(Deviance[,3], at = pos[3], add = T, 
        col = "lightblue", pch = 16, yaxt = "n", cex = 0.5)
legend("bottomleft",c("PD","FDP","DirMulti"),
       fill = c("darkred","darkblue","lightblue"), 
       cex = cex.legend, bty = "n", border = NA)
if(save_img)
  dev.off()




# GOF - multiple estimates ----------------------------------------------------------------
params_zipfs = list(0.9,1.02,2)
params_geom = list(0.85,0.9,0.95)
params_unif = list(NA)
params_negbin = list(c(1,0.003),c(5,0.003),c(1,0.01))
experiments = list("Zipfs" = params_zipfs,
                   "Geom" = params_geom,
                   "Uniform" = params_unif,
                   "NegBin" = params_negbin)
n = 500
M = 200 
M_max = 200
seed = 232131
set.seed(seed)
Mmax = M; Nrep = 50
AccCrv_length = 50
ii = 4; name = names(experiments)[[ii]]
jj = 1
params = experiments[[ii]][[jj]]
trim_params = sapply(params, get_first3digits, 4)
if(length(trim_params)>1)
  trim_params = paste0(trim_params[1],"_",trim_params[2])
ptrue = sim_generic_species(name,M,params)
ptrue = sort(ptrue, decreasing = TRUE)

Mplot = 100
par( mfrow = c(1,1), mar = c(3.5,4.25,1,1), mgp=c(2.75,1,0), bty = "l", las = 1, cex.lab = 1 )
plot(x = 1:Mplot, y = ptrue[1:Mplot], type = "p", pch = 16, cex = 0.5,
     main = paste0(name,", ",ii,", ",jj), xlab = "", ylab = "")
ptrue_mat = matrix(nrow = Nrep,ncol = M)
ptrue_mat = apply(ptrue_mat, 1, function(x) ptrue)

#a) Generate and fit multiple datasets
data_mat = SimData(n,ptrue_mat, seed) # n x Nrep
#b) Estimate parameters in each dataset
ParEst = t(apply(data_mat, 2, function(data){
  n_i = tabulate(data, nbins = M)
  idx_obs = which(n_i > 0)
  Kn = length(idx_obs)
  data_obs = n_i[idx_obs]
  
  res = matrix(nrow = 1, ncol = 5)
  colnames(res) = c("sigma-PD","theta-PD","Lambda-FDP","gamma-FDP","gamma-DirMulti")
  #a) PD
  # Param. estimation (PYP)
  start_params <- c(alpha = 0.5, theta = 1)
  fit <- optim(par = start_params, fn = llik_pyp, 
               n = n, Kn = Kn, data_obs = data_obs, # extra parameters
               method = "L-BFGS-B",
               lower = c(0, -1), upper = c(1-1e-10, Inf)) 
  res[1,1] = fit$par[1]
  res[1,2] = fit$par[2]
  
  #b) FDP
  # Param. estimation (FDP)
  start_params <- c(gamma = 0.1, Lambda = Kn)
  fit <- optim(par = start_params, fn = llik_FD, 
               n = n, Kn = Kn, data_obs = data_obs, M_max = M_max,# extra parameters
               method = "L-BFGS-B",
               lower = c(1e-5, 1e-5), upper = c(Inf, Inf)) 
  res[1,3] = fit$par[2]
  res[1,4] = fit$par[1]  
  
  #c) Dirichlet-Multinomial --> M = M
  # Param. estimation (DirMulti)
  start_params <- c(gamma = 0.5)
  fit <- optim(par = start_params, fn = llik_DirMult,
               n = n, M = M, data = n_i, # extra parameters
               method = "L-BFGS-B",
               lower = c(1e-10), upper = c(Inf))
  res[1,5] = fit$par[1]   
  
  res
})) # Nrep x 3 matrix
#c) Estimate probabilities for each dataset
ProbEst <- lapply(seq_len(nrow(ParEst)), function(i) {
  
  params <- ParEst[i, ]
  
  Mmax <- 3*M
  sigma    <- params[1]
  theta    <- params[2]
  Lambda   <- params[3]
  gamma    <- params[4]
  gamma_DM <- params[5]
  
  probPD  <- SimModel_generic(model = "PD",  Mmax = Mmax, Nrep = 1,
                              params = c(sigma, theta))
  
  probFDP <- SimModel_generic(model = "FDP", Mmax = Mmax, Nrep = 1,
                              params = c(gamma, Lambda))
  
  probDM  <- SimModel_generic(model = "DirMulti", Mmax = Mmax, Nrep = 1,
                              params = c(gamma_DM, M))
  
  probDM  <- rbind(probDM, matrix(0, nrow = Mmax - M, ncol = 1))
  
  cbind(probPD, probFDP, probDM)
})
ProbEst <- simplify2array(ProbEst)
ProbEst_qnt = apply(ProbEst, c(1,2), quantile, probs = c(0.025,0.5,0.975))
dim(ProbEst_qnt)
ProbEst_qnt[,1:5,]


## a - Envelop
xmax = M
# ymax = max( c(ptrue, max(gof_PD$Envelop_qnt), max(gof_FDP$Envelop_qnt), max(gof_DM$Envelop_qnt)) )
ymax = max( ptrue ) * 1.33
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos*1e3, 0)

img_name = paste0("img/GOF/",name,"_",trim_params,"_GOFEnv",".pdf")
if(save_img)
  pdf(img_name, width = width, height = height)
par( mfrow = c(1,3), mar = c(3.5,6,1,1), mgp=c(4.4,1,0), bty = "l", las = 1, cex.lab = cex.lab )

# PD
plot(0,0,main = " ", ylab = "Prob * 1000 (PD)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n",
     xlab = "",cex.axis = cex.axis)
axis( side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis) 
polygon( c(1:xmax, rev(1:xmax)),
         c(ProbEst_qnt[1,1:xmax,1], rev(ProbEst_qnt[3,1:xmax,1])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:xmax, y = ptrue[1:xmax], col = "darkred", pch = 16)

# FDP
plot(0,0,main = " ", ylab = "Prob * 1000 (FDP)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n",xlab = "",cex.axis = cex.axis)
axis( side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis) 
polygon( c(1:xmax, rev(1:xmax)),
         c(ProbEst_qnt[1,1:xmax,2], rev(ProbEst_qnt[3,1:xmax,2])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:xmax, y = ptrue[1:xmax], col = "darkred", pch = 16)

# DM
plot(0,0,main = " ", ylab = "Prob * 1000 (DiriMulti)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n",xlab = "",cex.axis = cex.axis)
axis( side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis) 
polygon( c(1:xmax, rev(1:xmax)),
         c(ProbEst_qnt[1,1:xmax,3], rev(ProbEst_qnt[3,1:xmax,3])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:xmax, y = ptrue[1:xmax], col = "darkred", pch = 16)


if(save_img)
  dev.off()

## b - Kn
Kn_true = apply(data_mat, 2, function(data){
  n_i = tabulate(data, nbins = M)
  idx_obs = which(n_i > 0)
  length(idx_obs)
})
Kn_est = apply(ProbEst, 2, function(pj_all){
  data_mat = SimData(n,pj_all) # n x Nrep
  apply(data_mat, 2, function(data){
    n_i = tabulate(data, nbins = M)
    idx_obs = which(n_i > 0)
    length(idx_obs)
  })
})
Kn_all = cbind(Kn_true,Kn_est)
colnames(Kn_all) = c("True","PD","FDP","DirMulti")


pos = c(1,3,5,7)
shift = c(0)
xpos = pos; xlabs = colnames(Kn_all)
ymax = max(Kn_all);ymin = min(Kn_all)
ylim_plot = c(ymin,ymax)

barplot_pos = c()
for(l in 1:length(pos)){
  barplot_pos = c(barplot_pos, pos[l]+shift)
}

img_name = paste0("img/GOF/",name,"_",trim_params,"_GOFKn",".pdf")
if(save_img)
  pdf(img_name, width = width, height = height)
par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(0,8), ylim = ylim_plot,
     xlab = "r", ylab = paste0("Kn"), 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
grid(lty = 1,lwd = 1, col = "gray90" )
boxplot(Kn_all[,1], at = pos[1], add = T, 
        col = "black", pch = 16, yaxt = "n", cex = 0.5)
boxplot(Kn_all[,2], at = pos[2], add = T, 
        col = "darkred", pch = 16, yaxt = "n", cex = 0.5)
boxplot(Kn_all[,3], at = pos[3], add = T, 
        col = "darkblue", pch = 16, yaxt = "n", cex = 0.5)
boxplot(Kn_all[,4], at = pos[4], add = T, 
        col = "lightblue", pch = 16, yaxt = "n", cex = 0.5)
legend("bottomleft",c("True","PD","FDP","DirMulti"),
       fill = c("black","darkred","darkblue","lightblue"), 
       cex = cex.legend, bty = "n", border = NA)
if(save_img)
  dev.off()

## c - Kn and Kn_r
rmax = 5
Knr_true = apply(data_mat, 2, function(data){
  n_i = tabulate(data, nbins = M)
  res = matrix(0,nrow = 1, ncol = rmax+1)
  colnames(res) = c("Kn",sapply(1:rmax, function(r) paste0("Kn_",r)))
  res[1,] = sapply(0:rmax, function(r){
    if(r == 0){
      idx_obs = which(n_i > 0)
    } else{
      idx_obs = which(n_i == r)
    }
    length(idx_obs)
  })
})
rownames(Knr_true) = c("Kn",sapply(1:rmax, function(r) paste0("Kn_",r))) # rmax+1 x Nrep

Nestimators = dim(ProbEst)[2]
Knr_list = vector("list", Nestimators)
names(Knr_list) = c("PD","FDP","DirMulti")
Knr_list = lapply(1:Nestimators, function(i){
  pj_all = ProbEst[,i,]
  data_mat = SimData(n,pj_all) # n x Nrep
  apply(data_mat, 2, function(data){
    n_i = tabulate(data, nbins = M)
    res = matrix(0,nrow = 1, ncol = rmax+1)
    colnames(res) = c("Kn",sapply(1:rmax, function(r) paste0("Kn_",r)))
    res[1,] = sapply(0:rmax, function(r){
      if(r == 0){
        idx_obs = which(n_i > 0)
      } else{
        idx_obs = which(n_i == r)
      }
      length(idx_obs)
    })
  })
})
Knr_list = lapply(Knr_list, function(x) {
  rownames(x) = c("Kn",sapply(1:rmax, function(r) paste0("Kn_",r))); x 
})

pos = c(1,3,5,7)
xpos = pos; xlabs = c("True","PD","FDP","DirMulti")

img_name = paste0("img/GOF/",name,"_",trim_params,"_GOF_Knr",".pdf")
if(save_img)
  pdf(img_name, width = width+6, height = height)
par( mfrow = c(2,3), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
for(ii in 1:(rmax+1)){
  Kn_plot = do.call(rbind, lapply(Knr_list, function(x) x[ii,])) 
  Kn_plot = rbind(Knr_true[ii,],Kn_plot)
  rownames(Kn_plot) = xlabs
  ymax = max(Kn_plot);ymin = min(Kn_plot)
  ylim_plot = c(ymin,ymax)
  ylab = ifelse(ii == 1, "Kn", paste0("Kn_",ii-1))
  

  plot(0,0,type = "n",xlim = c(0,8), ylim = ylim_plot,
       xlab = "r", ylab = ylab, 
       xaxt = "n", yaxt = "n",
       main = "",
       cex.lab = cex.lab, cex.axis = cex.axis)
  axis(side = 2, las = 1, cex.axis = cex.axis)
  axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
  grid(lty = 1,lwd = 1, col = "gray90" )
  boxplot(Kn_plot[1,], at = pos[1], add = T, 
          col = "black", pch = 16, yaxt = "n", cex = 0.5)
  boxplot(Kn_plot[2,], at = pos[2], add = T, 
          col = "darkred", pch = 16, yaxt = "n", cex = 0.5)
  boxplot(Kn_plot[3,], at = pos[3], add = T, 
          col = "darkblue", pch = 16, yaxt = "n", cex = 0.5)
  boxplot(Kn_plot[4,], at = pos[4], add = T, 
          col = "lightblue", pch = 16, yaxt = "n", cex = 0.5)
  # legend("bottomleft",c("PD","FDP","DirMulti"),
  #        fill = c("darkred","darkblue","lightblue"), 
  #        cex = cex.legend, bty = "n", border = NA)
}
if(save_img)
  dev.off()


## FDP - M and Mstar -------------------------------------------------------
colnames(ParEst) = c("sigma-PD","theta-PD","Lambda-FDP","gamma-FDP","gamma-DirMulti")
FDP_Lambda = ParEst[,3]
FDP_M = sapply(ParEst[,3], function(Lambda) rpois(n=1,lambda=Lambda))
FDP_Mstar = apply( cbind(ParEst,Knr_true[1,]), 1, function(params){
  Lambda = params[3]; gamma = params[4]; Kn = params[6]
  M_max = 500
  logExpMstar = compute_logV( Kn+1, n, gamma, Lambda, M_max ) - compute_logV( Kn, n, gamma, Lambda, M_max )
  ExpMstar = exp(logExpMstar) # Expected value
  ExpMstar
})
FDP_KnMstar = FDP_Mstar + Knr_true[1,]

pos = c(1,2)
FDP_M_plot = cbind(FDP_M,FDP_KnMstar)
colnames(FDP_M_plot) = c("M","Kn+Mstar")
ymax = max(FDP_M_plot);ymin = min(FDP_M_plot)
ylim_plot = c(ymin,ymax)
xpos = pos; xlabs = colnames(FDP_M_plot)

img_name = paste0("img/GOF/",name,"_",trim_params,"_GOF_FDP_M_Mstar",".pdf")
if(save_img)
  pdf(img_name, width = width, height = height)
par( mfrow = c(1,1), mar = c(3.5,6.5,1,1), mgp=c(5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
plot(0,0,type = "n",xlim = c(0.5,2.5), ylim = ylim_plot,
     xlab = "r", ylab = "(FDP)", 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
grid(lty = 1,lwd = 1, col = "gray90" )
boxplot(FDP_M_plot[,1], at = pos[1], add = T, 
        col = "steelblue", pch = 16, yaxt = "n", cex = 0.5)
boxplot(FDP_M_plot[,2], at = pos[2], add = T, 
        col = "#76EE00", pch = 16, yaxt = "n", cex = 0.5)
if(save_img)
  dev.off()







