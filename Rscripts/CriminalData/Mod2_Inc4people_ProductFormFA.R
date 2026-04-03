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
source("../../R/PFFAfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")


# Plot options ------------------------------------------------------------------
save_img = FALSE
width = 12; height = 6
cex.labels = 2
xmax = 140
# Load --------------------------------------------------------------------

load("RawDataInc.Rdat")
n = ncol(A)
Kn = nrow(A)
N_j = rowSums(A)
names(N_j) = as.character(1:Kn)
Nj_ordered = sort(N_j, decreasing = TRUE)

seed = 34231
set.seed(seed)

data = t(A)

eb_init_BB <- list(alpha = -1, s = 100, Nhat_prime = 50)
eb_known_BB <- list()
eb_params_obj_BB <- eb_params(model = "BB", init = eb_init_BB, known = eb_known_BB )

res = GibbsFA_eb(feature_matrix = data,
                 model = "NegBinBB_eb", type = "EFPF",
                 eb_params =  eb_params_obj_BB, 
                 var_fct = 1000)
res$mu0
res$fun_value
res$eb_params
res$alpha
res$theta

a_mle = res$alpha+1
b_mle = res$theta - a_mle
r_nb = res$n0
q_nb = 1 - res$mu0/res$var_fct

kappa_n = exp( lgamma(b_mle+n) + lgamma(a_mle+b_mle) - lgamma(b_mle) - lgamma(a_mle+b_mle+n) )
post_size = r_nb + Kn
post_Prfail = q_nb*kappa_n
ExpMstar = post_size * (post_Prfail)/(1 - post_Prfail)
ExpMstar

Mstar_mc = rnbinom(n = 10000,size = post_size, prob = 1 - post_Prfail)

pMstar = rep(0,max(Mstar_mc))
for(k in 1:(length(pMstar)-1)){
  pMstar[k] = dnbinom(x = (k-1), size = post_size, prob = 1 - post_Prfail)
}

bp1 = barplot(pMstar)

ypos = seq(0, max(pMstar), length.out = 3)
ynames = round(ypos*1000)


if(save_img)
  pdf("img/Mod3_Inc4people_MixBin_Mstar.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(1.5,0.5,0), mar = c(1,2.5,1,0), cex = 2)
barplot( height = pMstar, 
         names.arg = "", las = 2, col = "darkred", border = NA,
         main = " ", ylab = "Prob. * 1000", yaxt = "n" )
axis( side = 2, at = ypos, labels = ynames ,las = 1)
text( x = bp1, y = -0.001, 
      labels = 0:(length(pMstar)-1), 
      srt = 0, adj = 0.5, xpd = TRUE, cex = 1)
if(save_img)
  dev.off()

