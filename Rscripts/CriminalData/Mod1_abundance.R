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
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))

source("../../R/Rfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")

# Plot and colors ------------------------------------------------------------------

save_img = FALSE
width = 12; height = 6


# Load --------------------------------------------------------------------

load("Locale.Rdat")
Locale_temp = Locale[-which(Locale == "OUT" | Locale == "MISS")]
Nj_locale = table(Locale_temp)
Nj_locale = sort(Nj_locale, decreasing = TRUE)
bp1 = barplot(Nj_locale)

n = sum(Nj_locale)
Kn = length(Nj_locale)

seed = 34231
set.seed(seed)

# Check for uniformity ----------------------------------------------------

pval <- tryCatch( MultinomialTest(Nj = Nj_locale, M = NULL), error = function(e) 1 )
pval # small pvalue -> not uniform -> EB

# Upper bounds ------------------------------------------------------------
alfa <- alpha <- 0.05
Rmax = 100
var_prior <- 10

# (n_i,Mmax,M) are not known
Mguess = 20 # from juridical acts
n_i_guess = c(Nj_locale,rep(0,Mguess-Kn))
UB_all = UB_fit(n=n, Kn=Kn, n_i=n_i_guess,
                data_obs=Nj_locale,Mmax=NA,M=Mguess,
                var_prior=var_prior,Rmax=Rmax,alpha=alpha, 
                useMAP=TRUE,M_max=500,seed=seed)

## Freq. ------------------------------------------------------------
# ubfreq = ub_pain(n = n, Rmax = Rmax, alfa = alfa)
ubfreq = UB_all[1,2]
ubfreq

## PYP ------------------------------------------------------------
# # Param. estimation (PYP)
# start_params <- c(alpha = 0.1, theta = 1)
# fit <- optim(par = start_params, fn = llik_pyp, 
#              n = n, Kn = Kn, data_obs = Nj_locale, # extra parameters
#              method = "L-BFGS-B",
#              lower = c(0, -1), upper = c(1-1e-10, Inf)) 
# alpha_mle = fit$par[1]
# theta_mle = fit$par[2]

sigma_mle = UB_all[1,6]; theta_mle = UB_all[1,7]
c(sigma_mle,theta_mle)
# Plot
Nrep = 100; Natoms = 500
sim_PYP = r_SB(Nrep,Natoms,sigma_mle,theta_mle,seed)
sim_PYP = apply(sim_PYP, 1, sort, decreasing = TRUE)
pyp_qnt = apply(sim_PYP, 1, quantile, probs = c(0.025,0.5,0.975))

xmax = 20
ymax = 0.55#max(max(Nj_locale)/n , max(pyp_qnt) )
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)


if(save_img)
  pdf("img/Mod1_PD_paramfit.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(2.25,0.5,0), mar = c(2,3.25,1,1), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "Prob. (PD)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(pyp_qnt[1,1:xmax], rev(pyp_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:Kn, y = Nj_locale/n, col = "darkred", pch = 16, cex = 0.5)
segments(x0 = (1:Kn), x1 = (1:Kn), 
         y0 = rep(0,Kn), y1 = Nj_locale/n, col = "darkred", lwd = 3)
points(x = (1:xmax)+0.25, y = pyp_qnt[2,1:xmax], 
       col = "black", pch = 4, cex = 0.5)
segments(x0 = (1:xmax)+0.25, x1 = (1:xmax)+0.25,
         y0 = rep(0,xmax), y1 = pyp_qnt[2,1:xmax], 
         col = "black", lwd = 1, lty = 2)
if(save_img)
  dev.off()



# Upper bound (PYP)
# ubpyp = exp(compute_log_UBMarkov( Rmax, alpha_mle, theta_mle, Kn, n, alfa ))
ubpyp = UB_all[1,3]
ubpyp

## FDP ------------------------------------------------------------
M_max = 200
# Param. estimation (FD)
# start_params <- c(gamma = 0.1, Lambda = Kn)
# fit <- optim(par = start_params, fn = llik_FD, 
#              n = n, Kn = Kn, data_obs = Nj_locale, M_max = M_max,# extra parameters
#              method = "L-BFGS-B",
#              lower = c(1e-5, 1e-5), upper = c(Inf, Inf)) 
# gamma_mle = fit$par[1]
# Lambda_mle = fit$par[2]
gamma_FDP = UB_all[1,8]; Lambda_FDP = UB_all[1,9]
c(Lambda_FDP,gamma_FDP)

# Plot
Mub = 100
M_mc = rpois(n=Nrep,lambda = Lambda_FDP) + 1
sim_FD = matrix(0,nrow = Nrep, ncol = Mub)
ii = 1
for(ii in 1:Nrep){
  w = rgamma(n = M_mc[ii], shape = gamma_FDP, rate = 1); w = w/sum(w); w = sort(w, decreasing = TRUE)
  len = min(Mub,M_mc[ii])
  sim_FD[ii,1:len] = w[1:len]
}
sim_FD = apply(sim_FD, 1, sort, decreasing = TRUE)
FD_qnt = apply(sim_FD, 1, quantile, probs = c(0.025,0.5,0.975))

xmax = 20
ymax = 0.55 #max(max(Nj_locale)/n , max(FD_qnt) )
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

if(save_img)
  pdf("img/Mod1_FDP_paramfit.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(2.25,0.5,0), mar = c(2,3.25,1,1), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "Prob. (FDP)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(FD_qnt[1,1:xmax], rev(FD_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:Kn, y = Nj_locale/n, col = "darkred", pch = 16, cex = 0.5)
segments(x0 = (1:Kn), x1 = (1:Kn), 
         y0 = rep(0,Kn), y1 = Nj_locale/n, col = "darkred", lwd = 3)
points(x = (1:xmax)+0.25, y = FD_qnt[2,1:xmax], 
       col = "black", pch = 4, cex = 0.5)
segments(x0 = (1:xmax)+0.25, x1 = (1:xmax)+0.25,
         y0 = rep(0,xmax), y1 = FD_qnt[2,1:xmax], 
         col = "black", lwd = 1, lty = 2)
if(save_img)
  dev.off()

# Upper bound (FD)
# RmaxFD = 100
# ubFD = exp(compute_log_UBMarkov_FD( RmaxFD, gamma_mle, Lambda_mle, Kn, n, alfa, M_max ))
ubFD = UB_all[1,4]
ubFD

## DirMulti ------------------------------------------------------------
gamma_DM = UB_all[1,10]; gamma_DM

# Plot
sim_DM = matrix(0,nrow = Nrep, ncol = Mguess)
ii = 1
for(ii in 1:Nrep){
  w = rgamma(n = Mguess, shape = gamma_DM, rate = 1); w = w/sum(w); w = sort(w, decreasing = TRUE)
  sim_DM[ii,1:Mguess] = w[1:Mguess]
}
sim_DM = apply(sim_DM, 1, sort, decreasing = TRUE)
DM_qnt = apply(sim_DM, 1, quantile, probs = c(0.025,0.5,0.975))

xmax = 20
ymax = 0.55 #max(max(Nj_locale)/n , max(DM_qnt) )
ypos = seq(0, ymax, length.out = 5)
ylabs = round(ypos, 2)

if(save_img)
  pdf("img/Mod1_DirMulti_paramfit.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(2.25,0.5,0), mar = c(2,3.25,1,1), bty = "l", cex = 2)
plot(0,0,main = " ", ylab = "Prob. (Dir)", yaxt = "n", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
axis( side = 2, at = ypos, labels = ylabs, las = 1) 
polygon( c(1:xmax, rev(1:xmax)),
         c(DM_qnt[1,1:xmax], rev(DM_qnt[3,1:xmax])),
         col = "#FDE333",
         border = NA) # plot in-sample bands
points(x = 1:Kn, y = Nj_locale/n, col = "darkred", pch = 16, cex = 0.5)
segments(x0 = (1:Kn), x1 = (1:Kn), 
         y0 = rep(0,Kn), y1 = Nj_locale/n, col = "darkred", lwd = 3)
points(x = (1:xmax)+0.25, y = DM_qnt[2,1:xmax], 
       col = "black", pch = 4, cex = 0.5)
segments(x0 = (1:xmax)+0.25, x1 = (1:xmax)+0.25,
         y0 = rep(0,xmax), y1 = DM_qnt[2,1:xmax], 
         col = "black", lwd = 1, lty = 2)
if(save_img)
  dev.off()

# Upper bound (DirMulti)
ubDM = UB_all[1,5]
ubDM



# Missing mass analysis ------------------------------------------------------------------

## Freq ------------------------------------------------------------
N1 = length(which(Nj_locale == 1))
GT = N1/n
GT

Chao = SpadeR::ChaoSpecies(Nj_locale,datatype = "abundance")
coverage = SpadeR:::Chat.Ind(Nj_locale, m = sum(Nj_locale))   
coverage
Chao$Species_table

## PYP ------------------------------------------------------------
ExpMissMass_PYP = (theta_mle + sigma_mle*Kn)/(n + theta_mle)
ExpMissMass_PYP

xmax = 0.2;
xgrid = seq(0,xmax,length.out = 1000)
dMissMass_PYP = dbeta(x = xgrid, 
                      shape1 = theta_mle + sigma_mle*Kn,
                      shape2 = n - sigma_mle*Kn)

ymax = max(dMissMass_PYP) * 1.05
par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2.5,3.5,1,0), bty = "l")
plot(0,0,main = " ", ylab = "Missing mass (density)", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
points(x = xgrid, y = dMissMass_PYP,type = "l",lwd = 3, col = "black")
abline(v = ExpMissMass_PYP, col = "red", lty = 2, lwd = 2)
abline(v = qbeta(p = c(0.025,0.975), shape1 = theta_mle + sigma_mle*Kn, shape2 = n - sigma_mle*Kn ),
       col = "grey70", lty = 4, lwd = 2)



## FD ------------------------------------------------------------
Mstar_ub = 51

## Mstar
logExpMstar = compute_logV( Kn+1, n, gamma_FDP, Lambda_FDP, M_max ) - compute_logV( Kn, n, gamma_FDP, Lambda_FDP, M_max )
ExpMstar = exp(logExpMstar) # Expected value
ExpMstar

log_PMstar = rep(-Inf,Mstar_ub) # Whole distribution (log-scale)
log_PMstar = sapply(0:(Mstar_ub-1), function(m) log_qMpost(m,n,Kn,gamma_FDP,Lambda_FDP,M_max))

Mstar_MC = sample(0:(Mstar_ub-1), size = 10000, prob = exp(log_PMstar), replace = TRUE)
summary(Mstar_MC)
Mstar_tab = table(Mstar_MC)
bp1 = barplot(Mstar_tab)


if(save_img)
  pdf("img/Mod1_Mstar_barplot.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2.5,3.5,1,0), cex = 2)
barplot( height = Mstar_tab/length(Mstar_MC), 
         names.arg = "", las = 2, col = "darkred", border = NA,
         main = " ", ylab = "Mstar", yaxt = "n" )
axis( side = 2, at = seq(0, max(Mstar_tab/length(Mstar_MC)), by=0.01), las = 1)
text( x = bp1, y = -0.01, 
      labels = names(Mstar_tab), cex = 1,
      srt = 0, adj = 0.5, xpd = TRUE)
if(save_img)
  dev.off()

## Highest posterior density region
pMstar = Mstar_tab/length(Mstar_MC)
pMstar_ordered = sort(pMstar, decreasing = TRUE)
pMstar_ordered_csum = cumsum(pMstar_ordered)
if(save_img)
  pdf("img/Mod1_Mstar_HPD.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(2,0.5,0), mar = c(1.5,3,1,0), cex = 2)
barplot( height = pMstar_ordered_csum, 
         names.arg = "", las = 2, 
         col = c(rep("grey15",5),"grey40", c("grey70","grey70"), rep("grey90",4)),
         border = NA,
         main = " ", ylab = "HPD - Mstar", yaxt = "n" )
axis( side = 2, at = seq(0, 1, by=0.25), las = 1)
text( x = bp1, y = -0.05, 
      labels = names(pMstar_ordered_csum),
      srt = 0, adj = 0.5, xpd = TRUE)
abline( h = c(0.9,0.95,0.99), 
        col = c("grey15","grey40","grey70"), 
        lty = 2, lwd = 3 )
if(save_img)
  dev.off()
## Missing mass
log_ExpMissMass_vec = sapply(1:(Mstar_ub-1), function(m) {
  log(gamma_FDP) + log(m) - log(n + gamma_FDP*(Kn+m)) + log_PMstar[m+1]
})
ExpMissMass_FD = exp( log_stable_sum(log_ExpMissMass_vec, TRUE) ) # Expected value
ExpMissMass_FD

MissMass_FD = rbeta(n = length(Mstar_MC), 
                    shape1 = gamma_FDP*Mstar_MC,
                    shape2 = n + Kn*gamma_FDP)


dMissMass_FD = density(MissMass_FD)


ymax = max(dMissMass_FD$y) * 1.05
xmax = max(dMissMass_FD$x) * 1.05
par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2.5,3.5,1,0), bty = "l")
plot(0,0,main = " ", ylab = "Missing mass (density)", 
     xlim = c(0,xmax), ylim = c(0,ymax), type = "n")
points(x = dMissMass_FD$x, y = dMissMass_FD$y,
       type = "l",lwd = 3, col = "black")
abline(v = ExpMissMass_FD, col = "red", lty = 2, lwd = 2)
abline(v = quantile(MissMass_FD, probs = c(0.025,0.975)), 
       col = "grey70", lty = 4, lwd = 2)

