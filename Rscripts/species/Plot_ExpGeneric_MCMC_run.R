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
width = 12; height = 8
cex.labels <- cex.lab <- 1.5
cex.axis <- 1.5
cex.legend <- 1.5
mycol = c("darkorange","darkred","darkblue","lightblue")
mycol2 = c("black","lightblue")
lgd_names = c("Oracle","Freq","PD","FDP","Dir-Multi")

# Options -----------------------------------------------------------------
params_zipfs = list(0.9,1.02,2)
params_geom = list(0.85,0.9,0.95)
params_unif = list(NA)
params_negbin = list(c(1,0.003)) # list(c(1,0.003),c(5,0.003),c(1,0.01))
experiments = list("Zipfs" = params_zipfs,
                   "Geom" = params_geom,
                   "Uniform" = params_unif,
                   "NegBin" = params_negbin)

variance_prior_vec = c(1,10,100,1000)
alpha <- alfa <- 0.05
n = 500
Rmax = 100; RmaxFD = 50
Mmin_grid = 50; Mmax_grid = 1000
Mgrid = seq(Mmin_grid,Mmax_grid,by = 50); LMgrid = length(Mgrid)
Nexp = length(Mgrid)
M_max = 200

seed0 = 42
set.seed(seed0)

save_name_base = paste0("save/Species_MCMC_II_") 
save_name_base_cov = paste0("save/Species_MCMC_II_Cov_") 
img_fld = paste0("img/Species_") 


# Coverages ---------------------------------------------------------------

save_cov = FALSE
igrid = c(4)
ii = 4
cov_mat_print = vector("list",length = length(variance_prior_vec)); counter = 1
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
    var_prior = 10
    for(var_prior in variance_prior_vec){
      # Load
      # filename = paste0(save_name_base,name,"_nfix_",trim_params,".Rdat")
      filename = paste0(save_name_base,name,"_v_",var_prior,"_nfix_",trim_params,".Rdat")
      load(filename)
      
      # Coverage
      cov_mat <- do.call(rbind, lapply(ExpRes_list, function(x) {
        n <- nrow(x)
        if (n == 0) return(rep(NA_real_, ncol(x) - 1))
        counts <- colSums(x[, c(2,3,4,5), drop = FALSE] >= x[, 1])
        counts / n
      }))
      colnames(cov_mat) <- colnames(ExpRes_list[[1]])[c(2,3,4,5)]
      rownames(cov_mat) <- as.character(Mgrid)
      cov_mat_print[[counter]] = cov_mat
      counter = counter + 1
      if( any( apply(cov_mat,2,function(x) length(which(x<(1-alpha))))) )
        cat("\n To check: ",name," ",ii,"-",jj,"\n")
      
      if(save_cov)
        save(cov_mat, file = paste0(save_name_base_cov,name,"_v_",var_prior,"_nfix_",trim_params,".Rdat") )
        # save(cov_mat, file = paste0(save_name_base_cov,name,"_nfix_",trim_params,".Rdat") )
    }
    
  }
}

cov_mat_print

# Coverage - Table paper -----------------------------------
Tables = lapply(cov_mat_print, function(x) x[seq(2,LMgrid,by = 2),])

mat = cbind(Tables[[1]][,1],Tables[[2]][,1],Tables[[3]][,1],Tables[[4]][,1],
            Tables[[1]][,2],Tables[[2]][,2],Tables[[3]][,2],Tables[[4]][,2],
            Tables[[1]][,3],Tables[[2]][,3],Tables[[3]][,3],Tables[[4]][,3],
            Tables[[1]][,4],Tables[[2]][,4],Tables[[3]][,4],Tables[[4]][,4])
colnames(mat) = c(rep("Freq",4),rep("PD",4),rep("FDP",4),rep("DirMulti",4))
# Par Est --------------------------------------------------
ii = 4
igrid = c(4)
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
    
    PostMeans = vector("list",length = length(variance_prior_vec)); counter =  1
    for(var_prior in variance_prior_vec){
      # Load
      # filename = paste0(save_name_base,name,"_nfix_",trim_params,".Rdat")
      filename = paste0(save_name_base,name,"_v_",var_prior,"_nfix_",trim_params,".Rdat")
      load(filename)
      
      # Plot
      ParEst_mcmc = lapply(ExpRes_list, function(x) x[,c(6:10)] )
      ParEst_mcmc <- simplify2array(ParEst_mcmc) # Nrep x Nparams x LMgrid
      PostMeans[[counter]] = t(apply(ParEst_mcmc, c(2,3), mean)); counter = counter +1
      params_names = colnames(ExpRes_list[[1]])[6:10]
      for(ij in 1:dim(ParEst_mcmc)[2]){
        xx = ParEst_mcmc[,ij,]
        ## axis labels
        ymax = (11/10) * max(xx); ymin = (10/11) * min(xx)
        ylim_plot = c(ymin,ymax)
        ypos = seq(ymin,ymax,length.out = 5)
        ylabs = as.character(round(ypos,1))
        xmax = LMgrid+1; xmin = 0
        xlim_plot = c(0,xmax)
        xpos = 1:LMgrid
        xlabs = as.character(Mgrid)
        
        
        img_name = paste0(img_fld,name,"_MCMC","_v",var_prior,"_nfix_",trim_params,"_",params_names[ij],".pdf")
        if(save_img)
          pdf(img_name, width = width, height = height)
        par( mfrow = c(1,1), mar = c(3.5,6,1,1), mgp=c(4.5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
        plot(0,0,  yaxt = "n", xaxt = "n",
             xlab = "", ylab = paste0(var_prior," -- ",params_names[ij]),
             xlim = xlim_plot , ylim = ylim_plot, 
             main = paste0(" "),
             type = "n")
        grid(lty = 1,lwd = 1, col = "gray90" )
        axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
        axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
        mtext("M", side = 1, line = 2.5, cex = cex.axis)
        for(ii in 1:(LMgrid)){
          boxplot(xx[,ii], at = xpos[ii], add = T, 
                  col = "grey90", pch = 16, yaxt = "n", cex = 0.5)
        }
        if(save_img)
          dev.off()
        
      }
    }
    
  }
}


PostEst_gammaDM = do.call(cbind,lapply(PostMeans, function(x) x[seq(2,LMgrid,by = 2),5]))
rownames(PostEst_gammaDM) = Mgrid[seq(2,LMgrid,by = 2)]
colnames(PostEst_gammaDM) = variance_prior_vec
PostEst_gammaDM

# CI length ------------------------------------------------

ii = 4
igrid = c(4)
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
    
    for(var_prior in variance_prior_vec){
      # Load
      # filename = paste0(save_name_base,name,"_nfix_",trim_params,".Rdat")
      filename = paste0(save_name_base,name,"_v_",var_prior,"_nfix_",trim_params,".Rdat")
      load(filename)
      
      # Plot
      oracle = sapply(ExpRes_list, function(x) quantile(x[,1], 1-alpha))
      ExpRes_qnt = lapply(ExpRes_list, function(x) apply(x[,c(2:5)],2,quantile,probs = c(0.025,0.5,0.975),na.rm = TRUE) )
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
      
      
      img_name = paste0(img_fld,name,"_MCMC_","_v_",var_prior,"_nfix_",trim_params,".pdf")
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
      legend("topleft",lgd_names,
             fill = c("black",mycol), 
             cex = cex.legend, bty = "n", border = NA)
      if(save_img)
        dev.off()
    }

  }
}





# Brutta ------------------------------------------------------------------

var_prior = 1
filename = paste0(save_name_base,name,"_v_",var_prior,"_nfix_",trim_params,".Rdat")
load(filename)
ExpRes_list[[4]][1:5,c(6:10)]


var_prior = 10
filename = paste0(save_name_base,name,"_v_",var_prior,"_nfix_",trim_params,".Rdat")
load(filename)
ExpRes_list[[4]][1:5,c(6:10)]


var_prior = 1000
filename = paste0(save_name_base,name,"_v_",var_prior,"_nfix_",trim_params,".Rdat")
load(filename)
ExpRes_list[[4]][1:5,c(6:10)]






