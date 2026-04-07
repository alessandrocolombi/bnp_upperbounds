# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/features")
setwd(wd)

# Functions ---------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))

source("../../R/Rfunctions.R")
source("../../R/PFFAfunctions.R")
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")

# From BinomialCIs
source("../../../BinomialCIs/R/Rfunctions.R")
Rcpp::sourceCpp("../../../BinomialCIs/src/RcppFunctions.cpp")


# Custom functions --------------------------------------------------------

# Plot options ------------------------------------------------------------------
save_img = TRUE
width = 12; height = 6
cex.labels <- cex.lab <- 2
cex.axis <- 2
cex.legend <- 1.5
mycol = c("darkgreen","darkorange","darkred","darkblue","lightblue")
mycol2 = c("black","lightblue")
lgd_names = c("oracle","Bounded","Unbounded","IBP","MBP","FB")

# Options -----------------------------------------------------------------
params_zipfs  = list(0.25,0.5,1.2)
params_geom   = list(0.005,0.1,0.25)
params_const  = list(2,20,1000)
experiments   = list("Zipfs"   = params_zipfs,
                     "Geom"    = params_geom,
                     "Constant" = params_const)

alpha <- alfa <- 0.05
num_cores = 33 # <---
Nrep = 100 # <---
Rmax = 100; RmaxFD = 50
seed0 = 42
set.seed(seed0)
var_prior = 10
var_fct = 100
parallel = TRUE

# n fix -------------------------------------------------------------------

n = 2000
Mmin_grid = 100; Mmax_grid = 10000
Mgrid = seq(Mmin_grid,Mmax_grid,by = 500); LMgrid = length(Mgrid)
Nexp = length(Mgrid)
seeds = sample(1:999999, size = Nexp)

save_name_base = paste0("save/SS_features_nfix_") 
save_name_base_cov = paste0("save/SS_features_nfix_Cov") 
img_fld = paste0("img/SS_features_nfix_") 

## Coverage ----------------------------------------------------------------

save_cov = TRUE

ii = 1
igrid = c(1:3)
for(ii in igrid){
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
    filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
    load(filename)
    
    # Coverage
    cov_mat <- do.call(rbind, lapply(ExpRes_list, function(x) {
      x = x[,1:6]
      Nrep <- nrow(x)
      if (Nrep == 0) return(rep(NA_real_, ncol(x) - 1))
      counts <- colSums(x[, -1, drop = FALSE] >= x[, 1])
      counts / Nrep
    }))
    colnames(cov_mat) <- colnames(ExpRes_list[[1]])[2:6]
    rownames(cov_mat) <- as.character(Mgrid)
    
    if( any( apply(cov_mat,2,function(x) length(which(x<(1-alpha))))) )
      cat("\n To check: ",name," ",ii,"-",jj,"\n")
    
    if(save_cov)
      save(cov_mat, file = paste0(save_name_base_cov,name,"_",trim_params,".Rdat") )
  }
}

## Coverage - Table paper --------------------------------------------------


Tables = vector("list",3); counter = 1
ii = 1
for(ii in igrid){
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
    filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
    load(filename)
    load(paste0(save_name_base_cov,name,"_",trim_params,".Rdat"))
    cov_mat_all = cbind(cov_mat_all,cov_mat)
  }
  ord = c()
  for(j in 1:5){
    ord = c(ord, seq(j,5*Ncases,by = 5))
  }
  cov_mat_all = cov_mat_all[,ord]
  cov_mat_all = round(cov_mat_all,3)
  Tables[[counter]] = cov_mat_all
  counter = counter + 1
}
idx_rows = seq(1,LMgrid,by = 3)
Tables = lapply(Tables, function(x) x[idx_rows,])
Tables[[1]]#[,-c(4,8,12,16)]
Tables[[2]]
Tables[[3]]

## CI length ---------------------------------------------------------------
ii = 1
igrid = c(1:3)
for(ii in igrid){
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
    filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
    load(filename)
    
    # Plot
    oracle = sapply(ExpRes_list, function(x) quantile(x[,1], 1-alpha))
    ExpRes_qnt = lapply(ExpRes_list, function(x) apply(x[,c(2:6)],2,quantile,probs = c(0.025,0.5,0.975)) )
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
    
    
    img_name = paste0(img_fld,name,"_",trim_params,".pdf")
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
      polygon( c(Mgrid, rev(Mgrid)),
               c(ExpRes_qnt[1,ij,], rev(ExpRes_qnt[3,ij,])),
               col = scales::alpha(mycol[ij], 0.25),
               border = NA) # plot in-sample bands
    }
    # legend("topright",lgd_names,
    #        fill = c("black",mycol), 
    #        cex = cex.legend, bty = "n", border = NA)
    if(save_img)
      dev.off()
  }
}

## Par Est --------------------------------------------------

# Aiutati --> non salvare tutte queste immagini
ii = 1
igrid = c(1:3)
for(ii in igrid){
  
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 2
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = sapply(params, get_first3digits, 4)
    if(length(trim_params)>1)
      trim_params = paste0(trim_params[1],"_",trim_params[2])
    
    # Load
    filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
    load(filename)
    
    # Plot
    ParEst_mcmc = lapply(ExpRes_list, function(x) x[,c(7:16)] ) # Kn in included here
    ParEst_mcmc <- simplify2array(ParEst_mcmc) # Nrep x Nparams x LMgrid
    ExpMstar = apply(ParEst_mcmc, 3, function(x){
      apply(x,1,function(xx){
        a = xx[4]; b = xx[5]; m = xx[6]; q = xx[7]; Kn = xx[10]
        kappa_n = exp( lgamma(b+n)-lgamma(b)+lgamma(a+b)-lgamma(a+b+n) )
        (m+Kn) * q*kappa_n/(1 - 1*kappa_n)
      })
    }) # Nrep x LMgrid
    # ParEst_mcmc <- array(c(ParEst_mcmc, ExpMstar),dim = c(Nrep, 7, LMgrid))
    params_names = c(colnames(ExpRes_list[[1]])[7:16],"Kn+Mstar")
    for(ij in 1:(length(params_names))){
      if(ij < length(params_names) ){
        xx = ParEst_mcmc[,ij,]
      } else if(ij == length(params_names)){
        xx = ExpMstar
      } else{
        stop("non va")
      }
      ## axis labels
      ymax = (11/10) * max(xx); ymin = (10/11) * min(xx)
      ylim_plot = c(ymin,ymax)
      ypos = seq(ymin,ymax,length.out = 5)
      ylabs = as.character(round(ypos,1))
      xmax = LMgrid+1; xmin = 0
      xlim_plot = c(0,xmax)
      xpos = 1:LMgrid
      xlabs = as.character(Mgrid)
      
      
      img_name = paste0(img_fld,name,"_",trim_params,"_",params_names[ij],".pdf")
      if(save_img)
        pdf(img_name, width = width, height = height)
      par( mfrow = c(1,1), mar = c(3.5,6,1,1), mgp=c(4.5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
      plot(0,0,  yaxt = "n", xaxt = "n",
           xlab = "", ylab = paste0(params_names[ij]),
           xlim = xlim_plot , ylim = ylim_plot, 
           main = paste0(" "),
           type = "n")
      grid(lty = 1,lwd = 1, col = "gray90" )
      axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
      axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
      mtext("M", side = 1, line = 2.5, cex = cex.axis)
      for(iii in 1:(LMgrid)){
        boxplot(xx[,iii], at = xpos[iii], add = T, 
                col = "grey90", pch = 16, yaxt = "n", cex = 0.5)
      }
      if(save_img)
        dev.off()
      
    }
  }
}



# M fix -------------------------------------------------------------------

M = 5000
Nmin_grid = 500; Nmax_grid = 10000
Ngrid = seq(Nmin_grid,Nmax_grid,by = 500); LNgrid = length(Ngrid)
Nexp = length(Ngrid)
seeds = sample(1:999999, size = Nexp)

save_name_base = paste0("save/SS_features_Mfix_") 
save_name_base_cov = paste0("save/SS_features_Mfix_Cov_") 
img_fld = paste0("img/SS_features_Mfix_") 



## Coverage ----------------------------------------------------------------

save_cov = TRUE

ii = 1
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
    filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
    load(filename)
    
    # Coverage
    cov_mat <- do.call(rbind, lapply(ExpRes_list, function(x) {
      x = x[,1:6]
      Nrep <- nrow(x)
      if (Nrep == 0) return(rep(NA_real_, ncol(x) - 1))
      counts <- colSums(x[, -1, drop = FALSE] >= x[, 1])
      counts / Nrep
    }))
    colnames(cov_mat) <- colnames(ExpRes_list[[1]])[2:6]
    rownames(cov_mat) <- as.character(Ngrid)
    
    if( any( apply(cov_mat,2,function(x) length(which(x<(1-alpha))))) )
      cat("\n To check: ",name," ",ii,"-",jj,"\n")
    
    if(save_cov)
      save(cov_mat, file = paste0(save_name_base_cov,name,"_",trim_params,".Rdat") )
  }
}

## Coverage - Table paper --------------------------------------------------


Tables = vector("list",3); counter = 1
ii = 1
for(ii in 1:3){
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 2
  cov_mat_all = matrix(nrow = LNgrid, ncol = 0)
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = sapply(params, get_first3digits, 4)
    if(length(trim_params)>1)
      trim_params = paste0(trim_params[1],"_",trim_params[2])
    
    # Load
    filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
    load(filename)
    load(paste0(save_name_base_cov,name,"_",trim_params,".Rdat"))
    cov_mat_all = cbind(cov_mat_all,cov_mat)
  }
  ord = c()
  for(j in 1:5){
    ord = c(ord, seq(j,5*Ncases,by = 5))
  }
  cov_mat_all = cov_mat_all[,ord]
  cov_mat_all = round(cov_mat_all,3)
  Tables[[counter]] = cov_mat_all
  counter = counter + 1
}
idx_rows = c(1,2,3,5,7,9)
Tables = lapply(Tables, function(x) x[idx_rows,])
Tables[[1]]
Tables[[2]]
Tables[[3]]

## CI length ---------------------------------------------------------------
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
    filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
    load(filename)
    
    # Plot
    oracle = sapply(ExpRes_list, function(x) quantile(x[,1], 1-alpha))
    ExpRes_qnt = lapply(ExpRes_list, function(x) apply(x[,c(2:5)],2,quantile,probs = c(0.025,0.5,0.975)) )
    ExpRes_qnt <- simplify2array(ExpRes_qnt) # 3 x 5 x LNgrid
    
    ## axis labels
    ymax = (11/10) * max(ExpRes_qnt,oracle); ymin = (10/11) * min(ExpRes_qnt,oracle)
    # ymax = 17.5*1e-3; ymin = 11.5*1e-3
    ylim_plot = c(ymin,ymax)
    ypos = seq(ymin,ymax,length.out = 5)
    ylabs = as.character(round(ypos*1e3,0))
    xmax = max(Ngrid); xmin = min(Ngrid)
    xlim_plot = c(0,xmax)
    xpos = Ngrid
    xlabs = as.character(Ngrid)
    
    
    img_name = paste0(img_fld,name,"_",trim_params,".pdf")
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
    mtext("n", side = 1, line = 2.5, cex = cex.axis)
    points( x = Ngrid, y = oracle, 
            type = "l", lwd = 5, col = "black" )
    for(ij in 1:dim(ExpRes_qnt)[2]){
      points( x = Ngrid, y = ExpRes_qnt[2,ij,], 
              type = "l", lwd = 5, col = mycol[ij] )
      polygon( c(Ngrid, rev(Ngrid)),
               c(ExpRes_qnt[1,ij,], rev(ExpRes_qnt[3,ij,])),
               col = scales::alpha(mycol[ij], 0.25),
               border = NA) # plot in-sample bands
    }
    legend("topright",lgd_names,
           fill = c("black",mycol), 
           cex = cex.legend, bty = "n", border = NA)
    if(save_img)
      dev.off()
  }
}

## Par Est --------------------------------------------------

# Aiutati --> non salvare tutte queste immagini
ii = 1
igrid = c(1:4)
for(ii in igrid){
  
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
    filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
    load(filename)
    
    # Plot
    ParEst_mcmc = lapply(ExpRes_list, function(x) x[,c(6:11)] ) # Kn in included here
    ParEst_mcmc <- simplify2array(ParEst_mcmc) # Nrep x Nparams x LNgrid
    params_names = c(colnames(ExpRes_list[[1]])[6:11])
    for(ij in 1:6){
      if(ij < 7){
        xx = ParEst_mcmc[,ij,]
      } else if(ij == 7){
        xx = ExpMstar
      } else{
        stop("non va")
      }
      ## axis labels
      ymax = (11/10) * max(xx); ymin = (10/11) * min(xx)
      ylim_plot = c(ymin,ymax)
      ypos = seq(ymin,ymax,length.out = 5)
      ylabs = as.character(round(ypos,1))
      xmax = LNgrid+1; xmin = 0
      xlim_plot = c(0,xmax)
      xpos = 1:LNgrid
      xlabs = as.character(Ngrid)
      
      
      img_name = paste0(img_fld,name,"_",trim_params,"_",params_names[ij],".pdf")
      if(save_img)
        pdf(img_name, width = width, height = height)
      par( mfrow = c(1,1), mar = c(3.5,6,1,1), mgp=c(4.5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
      plot(0,0,  yaxt = "n", xaxt = "n",
           xlab = "", ylab = paste0(params_names[ij]),
           xlim = xlim_plot , ylim = ylim_plot, 
           main = paste0(" "),
           type = "n")
      grid(lty = 1,lwd = 1, col = "gray90" )
      axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
      axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
      mtext("M", side = 1, line = 2.5, cex = cex.axis)
      for(iii in 1:(LNgrid)){
        boxplot(xx[,iii], at = xpos[iii], add = T, 
                col = "grey90", pch = 16, yaxt = "n", cex = 0.5)
      }
      if(save_img)
        dev.off()
      
    }
  }
}





