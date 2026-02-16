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
params_geom = list(0.5,1-0.8,1-0.95)
params_unif = list(NA)
params_negbin = list(c(1,0.003),c(5,0.003),c(1,0.01))
experiments = list("Zipfs" = params_zipfs,
                   "Geom" = params_geom,
                   "Uniform" = params_unif,
                   "NegBin" = params_negbin)
Nrep = 200 # <---
n = 500
M = 400 
M_max = 200
seed = 232131
set.seed(seed)
Mmax = M; Nrep = 50
AccCrv_length = 50
ii = 2; name = names(experiments)[[ii]]
jj = 2
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

# Plots

## a - Envelop
xmax = 500
ymax = max( c(ptrue, max(gof_PD$Envelop_qnt), max(gof_FDP$Envelop_qnt)) )
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

img_name = paste0("img/GOF/",name,"_",trim_params,"_GOFEnv",".pdf")
if(save_img)
  pdf(img_name, width = width, height = height)
par( mfrow = c(1,2), mar = c(3.5,6,1,1), mgp=c(4.4,1,0), bty = "l", las = 1, cex.lab = cex.lab )

# PD
plot(0,0,main = " ", ylab = "#(PD)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n",
     xlab = "",cex.axis = cex.axis)
axis( side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis) 
polygon( c(1:xmax, rev(1:xmax)),
         c(gof_PD$Envelop_qnt[1,1:xmax], rev(gof_PD$Envelop_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:xmax, y = ptrue[1:xmax], col = "darkred", pch = 16)

# FDP
plot(0,0,main = " ", ylab = "#(FDP)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n",xlab = "",cex.axis = cex.axis)
axis( side = 2, at = ypos, labels = ylabs, las = 1, cex.axis = cex.axis) 
polygon( c(1:xmax, rev(1:xmax)),
         c(gof_FDP$Envelop_qnt[1,1:xmax], rev(gof_FDP$Envelop_qnt[3,1:xmax])),
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
     xlab = "n", ylab = paste0("Freq.Rare"), 
     xaxt = "n", yaxt = "n",
     main = "",
     cex.lab = cex.lab, cex.axis = cex.axis)
axis(side = 2, las = 1, cex.axis = cex.axis)
axis(side = 1, at = xpos, labels = xlabs, las = 1, cex.axis = cex.axis)
points(x = ngrid, y = gof_true$AccCrv[2,], type = "l", lwd = 5, col = "black")
points(x = ngrid, y = gof_PD$AccCrv[2,], type = "l", lwd = 5, col = "darkred")
points(x = ngrid, y = gof_FDP$AccCrv[2,], type = "l", lwd = 5, col = "darkblue")
if(save_img)
  dev.off()

