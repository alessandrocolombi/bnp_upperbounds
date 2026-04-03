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

# Plot options ------------------------------------------------------------------
save_img = FALSE
width = 12; height = 6
cex.labels <- cex.lab <- 2
cex.axis <- 2
cex.legend <- 1.5
mycol = c("darkblue","darkred","darkorange","lightblue")
mycol2 = c("black","lightblue")
lgd_names = c("Freq","PYP","FDP","Dir-Multi")

# Load --------------------------------------------------------------------

load("Locale.Rdat")
Locale_temp = Locale[-which(Locale == "OUT" | Locale == "MISS")]
data = Locale_temp
Nj_locale = table(Locale_temp)
Nj_locale = sort(Nj_locale, decreasing = TRUE)

n = sum(Nj_locale)
Kn = length(Nj_locale)


# Options  --------------------------------------------------------
eps_grid = seq(0.001,0.2,length.out =  34*5)
cov_grid = 1 - eps_grid
alpha <- alfa <- 0.05
M_max = 200
nstart = 10
var_prior = 10
Mguess = 20

seed0 = 4224
num_cores = 33 # <---
Nrep = 50


# Run) Mmax-based  --------------------------------------------------------
# res = SRabu_grid( eps_grid=eps_grid, data=data, nstart=nstart,
#                   Mguess=Mguess, var_prior=var_prior,
#                   Nrep=Nrep, num_cores=num_cores, seed0=seed0, alpha=alpha, M_max=M_max)
# save(res, file = "save/Mod1Abu_SRMmax.Rdat")

# Run) Coverage-based  --------------------------------------------------------
# res_cov = SRabu_cov_grid( cov_grid, data, nstart, Nrep, num_cores, seed0)
# save(res_cov, file = "save/Mod1Abu_SRcov.Rdat")

# Plot --------------------------------------------------------------------
stop_here = TRUE
ltype = c(1,1,1,1)
ygrids = vector("list",4)
ygrids[[1]]<-ygrids[[2]]<-ygrids[[3]]<-eps_grid
ygrids[[4]]<- eps_grid #(1-cov_grid)

if(!stop_here){
  load("save/Mod1Abu_SRMmax.Rdat")
  # load("save/Mod1Abu_SRcov.Rdat")
  
  res_list = lapply(res, function(x) apply(x,2,quantile,probs = c(0.025,0.5,0.975)))
  # res_cov_list = lapply(res_cov, function(x) apply(x,2,quantile,probs = c(0.025,0.5,0.975)))
  res_arr <- simplify2array(res_list) # 3 x 3 x length(eps_grid)
  # res_cov_arr <- simplify2array(res_cov_list) # 3 x 4 x length(eps_grid)
  # res_all <- array(
  #   do.call(cbind, lapply(seq_len(dim(res_arr)[3]), function(k)
  #     cbind(res_arr[,,k], res_cov_arr[,,k])
  #   )),
  #   dim = c(3, 4, 70)
  # )# 3 x 4 x length(eps_grid)
  res_all = res_arr
  
  
  
  ## Plot
  ## axis labels
  ymax = (11/10) * max(res_all); 
  ymin = 0 #(10/11) * min(res_all)
  ylim_plot = c(ymin,ymax)
  ypos = seq(ymin,ymax,length.out = 5)
  ylabs = as.character(round(ypos,0))
  xmax = max(eps_grid); xmin = min(eps_grid)
  xlim_plot = c(0,xmax)
  xpos = seq(xmin,xmax,length.out = 5)
  xlabs = as.character(round(xpos,2))
  
  if(save_img)
    pdf("img/Mod1Abu_StopR_grid.pdf", width = width, height = height)
  par( mfrow = c(1,1), mar = c(3.5,4.25,1,1), mgp=c(2.75,1,0), bty = "l", las = 1, cex.lab = cex.lab )
  plot(0,0,  yaxt = "n", xaxt = "n",
       xlab = "", ylab = "N stop",
       xlim = xlim_plot , ylim = ylim_plot, 
       main = paste0(" "),
       type = "n")
  grid(lty = 1,lwd = 1, col = "gray90" )
  axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
  axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
  mtext(expression(epsilon), side = 1, line = 2.5, cex = cex.axis)
  for(ij in 1:dim(res_all)[2]){
    points( x = eps_grid, y = res_all[2,ij,], 
            type = "l", lwd = 5, col = mycol[ij] )
    polygon( c(eps_grid, rev(eps_grid)),
             c(res_all[1,ij,], rev(res_all[3,ij,])),
             col = scales::alpha(mycol[ij], 0.25),
             border = NA) # plot in-sample bands
  }
  legend("topright",lgd_names,
         fill = c(mycol), 
         cex = cex.legend, bty = "n", border = NA)
  if(save_img)
    dev.off()
}

