# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
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
mycol = c("darkorange","darkred","darkblue","lightblue")
mycol2 = c("black","lightblue")
lgd_names = c("Oracle","Freq","PD","FDP","Dir-Multi")

# Options -----------------------------------------------------------------
params_zipfs = list(0.9,1.02,2)
params_geom = list(0.85,0.9,0.95)
params_unif = list(NA)
params_negbin = list(c(1,0.003),c(5,0.003),c(1,0.01))
experiments = list("Zipfs" = params_zipfs,
                   "Geom" = params_geom,
                   "Uniform" = params_unif,
                   "NegBin" = params_negbin)

alpha <- alfa <- 0.05
Nrep = 200 # <---
n = 500
Rmax = 100; RmaxFD = 50
Mmin_grid = 50; Mmax_grid = 1000
Mgrid = seq(Mmin_grid,Mmax_grid,by = 50); LMgrid = length(Mgrid)
Nexp = length(Mgrid)
M_max = 200

seed0 = 42
set.seed(seed0)

save_name_base = paste0("save/Species_") 
save_name_base_cov = paste0("save/Species_Cov_") 
img_fld = paste0("img/Species_") 


# Coverages ---------------------------------------------------------------

save_cov = FALSE

ii = 4
for(ii in 1:length(experiments)){
  
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 3
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = sapply(params, get_first3digits, 4)
    if(length(trim_params)>1)
      trim_params = paste0(trim_params[1],"_",trim_params[2])
    
    # Load
    filename = paste0(save_name_base,name,"_nfix_",trim_params,".Rdat")
    load(filename)
    
    # Coverage
    cov_mat <- do.call(rbind, lapply(ExpRes_list, function(x) {
      n <- nrow(x)
      if (n == 0) return(rep(NA_real_, ncol(x) - 1))
      counts <- colSums(x[, -1, drop = FALSE] >= x[, 1])
      counts / n
    }))
    colnames(cov_mat) <- colnames(ExpRes_list[[1]])[-1]
    rownames(cov_mat) <- as.character(Mgrid)
    
    if( any( apply(cov_mat,2,function(x) length(which(x<(1-alpha))))) )
      cat("\n To check: ",name," ",ii,"-",jj,"\n")
    
    if(save_cov)
      save(cov_mat, file = paste0(save_name_base_cov,name,"_nfix_",trim_params,".Rdat") )
  }
}


# Coverage - Table paper --------------------------------------------------

# Zipfs and Geometric
ii = 2
for(ii in 1:2){
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 2
  cov_mat_all = matrix(nrow = LMgrid, ncol = 0)
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = sapply(params, get_first3digits, 4)
    if(length(trim_params)>1)
      trim_params = paste0(trim_params[1],"_",trim_params[2])
    
    # Load
    filename = paste0(save_name_base,name,"_nfix_",trim_params,".Rdat")
    load(filename)
    load(paste0(save_name_base_cov,name,"_nfix_",trim_params,".Rdat"))
    cov_mat_all = cbind(cov_mat_all,cov_mat)
  }
  cov_mat_all = cov_mat_all[,c(1,5,9,2,6,10,3,7,11,4,8,12)]
  cov_mat_all = cov_mat_all[seq(1,20,by=3),]
  round(cov_mat_all,2)
}

# Uniform
ii = 3
for(ii in 3:3){
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 1
  cov_mat_all = matrix(nrow = LMgrid, ncol = 0)
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = sapply(params, get_first3digits, 4)
    if(length(trim_params)>1)
      trim_params = paste0(trim_params[1],"_",trim_params[2])
    
    # Load
    filename = paste0(save_name_base,name,"_nfix_",trim_params,".Rdat")
    load(filename)
    load(paste0(save_name_base_cov,name,"_nfix_",trim_params,".Rdat"))
    cov_mat_all = cbind(cov_mat_all,cov_mat)
  }
  cov_mat_all = cov_mat_all[seq(1,20,by=3),]
  round(cov_mat_all,3)
}

# Negative Binomial
ii = 4
for(ii in 4:4){
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 1
  cov_mat_all = matrix(nrow = LMgrid, ncol = 0)
  for(jj in c(1,3)){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = sapply(params, get_first3digits, 4)
    if(length(trim_params)>1)
      trim_params = paste0(trim_params[1],"_",trim_params[2])
    
    # Load
    filename = paste0(save_name_base,name,"_nfix_",trim_params,".Rdat")
    load(filename)
    load(paste0(save_name_base_cov,name,"_nfix_",trim_params,".Rdat"))
    cov_mat_all = cbind(cov_mat_all,cov_mat)
  }
  cov_mat_all = cov_mat_all[,c(5,1,6,2,7,3,8,4)]
  cov_mat_all = cov_mat_all[seq(1,20,by=3),]
  round(cov_mat_all,3)
}
# CI length ---------------------------------------------------------------

ii = 3
for(ii in 1:length(experiments)){
  
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 1
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = sapply(params, get_first3digits, 4)
    if(length(trim_params)>1)
      trim_params = paste0(trim_params[1],"_",trim_params[2])
    
    # Load
    filename = paste0(save_name_base,name,"_nfix_",trim_params,".Rdat")
    load(filename)
    
    # Plot
    oracle = sapply(ExpRes_list, function(x) quantile(x[,1], 1-alpha))
    ExpRes_qnt = lapply(ExpRes_list, function(x) apply(x[,c(2:5)],2,quantile,probs = c(0.025,0.5,0.975)) )
    ExpRes_qnt <- simplify2array(ExpRes_qnt) # 3 x 5 x LMgrid
    
    ## axis labels
    ymax = (11/10) * max(ExpRes_qnt,oracle); ymin = (10/11) * min(ExpRes_qnt,oracle)
    # ymax = 17.5*1e-3; ymin = 11.5*1e-3
    ylim_plot = c(ymin,ymax)
    ypos = seq(ymin,ymax,length.out = 5)
    ylabs = as.character(round(ypos*1e3,0))
    xmax = max(Mgrid); xmin = min(Mgrid)
    xlim_plot = c(0,xmax)
    xpos = Mgrid
    xlabs = as.character(Mgrid)
    
    
    img_name = paste0(img_fld,name,"_nfix_",trim_params,".pdf")
    if(save_img)
      pdf(img_name, width = width, height = height)
    par( mfrow = c(1,1), mar = c(3.5,4.25,1,1), mgp=c(2.75,1,0), bty = "l", las = 1, cex.lab = cex.lab )
    plot(0,0,  yaxt = "n", xaxt = "n",
         xlab = "", ylab = "1000 * bound",
         xlim = xlim_plot , ylim = ylim_plot, 
         main = paste0(" "),
         type = "n")
    grid(lty = 1,lwd = 1, col = "gray90" )
    axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
    axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
    mtext("M", side = 1, line = 2.5, cex = cex.axis)
    points( x = Mgrid, y = oracle, 
            type = "l", lwd = 5, col = "black" )
    for(ij in 1:dim(ExpRes_qnt)[2]){
      points( x = Mgrid, y = ExpRes_qnt[2,ij,], 
              type = "l", lwd = 5, col = mycol[ij] )
      # points( x = Mgrid, y = ExpRes_qnt[1,ij,],
      # type = "l", lwd = 3, lty = 2, col = mycol[ij] )
      # points( x = Mgrid, y = ExpRes_qnt[3,ij,], 
      #         type = "l", lwd = 3, lty = 2, col = mycol[ij] )
    }
    legend("topleft",lgd_names,
           fill = c("black",mycol), 
           cex = cex.legend, bty = "n", border = NA)
    if(save_img)
      dev.off()
  }
}












