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


compute_ab_beta = function(m, v, tol = sqrt(.Machine$double.eps) ) 
{
  if (length(m) != 1 || length(v) != 1) 
    stop("m and v must be scalars")
  if (!is.finite(m) || !is.finite(v)) 
    stop("m and v must be finite")
  if (m <= 0 || m >= 1) 
    stop("m must satisfy 0 < m < 1")
  if (v <= 0) 
    stop("v must satisfy v > 0")
  
  vmax = m * (1 - m)
  
  if (v >= vmax) 
    stop("Need v < m*(1-m) for a proper Beta distribution")
  
  # detect near-boundary regime: a+b very close to 0
  if ((vmax - v) <= tol * vmax) {
    warning("v is extremely close to m*(1-m): a and b are near 0 and numerically unstable")
  }
  
  kappa = (vmax - v) / v
  a = m * kappa
  b = (1 - m) * kappa
  
  c(a = a, b = b)
}


# Custom functions --------------------------------------------------------

# Plot options ------------------------------------------------------------------
save_img = FALSE
width = 12; height = 6
cex.labels <- cex.lab <- 2
cex.axis <- 2
cex.legend <- 1.5
mycol = c("darkgreen","darkorange","darkred","darkblue","lightblue")
mycol2 = c("black","lightblue")
lgd_names = c("oracle","Bounded","Unbounded","IBP","MBP","FB")

# Options -----------------------------------------------------------------
params_zipfs  = list(0.85,1.02,1.2)
params_geom   = list(0.005,0.1,0.25)
params_const  = list(2,1000,5000)
experiments   = list("Zipfs"   = params_zipfs,
                     "Geom"    = params_geom,
                     "Constant" = params_const)

alpha <- alfa <- 0.05
num_cores = 33 # <---
Nrep = 500 # <---
Rmax = 100; RmaxFD = 50
seed0 = 42
set.seed(seed0)
var_prior = 10
var_fct = 100
parallel = TRUE

format_cov_triplets = function(x, digits = 2) {
  apply(x, 1, function(row_vals) {
    vals = formatC(row_vals, format = "f", digits = digits)
    paste(vals, collapse = " / ")
  })
}

write_latex_cov_table = function(tab, row_values, row_name, caption = "", label, file) {
  nr = nrow(tab)
  nc = ncol(tab)
  n_methods = 5
  n_cases = nc / n_methods
  if (abs(n_cases - round(n_cases)) > .Machine$double.eps^0.5) {
    stop("Table does not have a multiple of 5 columns.")
  }
  n_cases = as.integer(n_cases)
  method_labels = c("Bdd", "Ubd", "IBP", "MBP", "FB")
  body_cols = lapply(seq_len(n_methods), function(j) {
    idx = seq(j, nc, by = n_methods)
    format_cov_triplets(tab[, idx, drop = FALSE])
  })
  body_df = data.frame(
    row_value = row_values,
    Bdd = body_cols[[1]],
    Ubd = body_cols[[2]],
    IBP = body_cols[[3]],
    MBP = body_cols[[4]],
    FB = body_cols[[5]],
    stringsAsFactors = FALSE
  )
  lines = c(
    "\\begin{table}[ht!]",
    "    \\centering",
    "    \\small",
    "    \\begin{tabular}{r c c c c c}",
    "    \\hline",
    paste0("    $", row_name, "$ & ", paste(method_labels, collapse = " & "), " \\\\"),
    "    \\hline"
  )
  for (i in seq_len(nr)) {
    lines = c(lines, paste0(
      "    ", body_df$row_value[i], " & ",
      paste(body_df[i, method_labels], collapse = " & "),
      " \\\\"
    ))
  }
  lines = c(
    lines,
    "    \\hline",
    "    \\end{tabular}",
    paste0("    \\caption{", caption, "}"),
    paste0("    \\label{", label, "}"),
    "\\end{table}"
  )
  writeLines(lines, con = file)
  cat(paste(lines, collapse = "\n"), "\n\n")
}

make_param_tag = function(params) {
  if (length(params) == 1 && is.na(params)) {
    return("NA")
  }
  fmt_one = function(x) {
    sx = format(x, scientific = FALSE, trim = TRUE)
    sx = gsub("\\.", "p", sx)
    sx = gsub("-", "m", sx)
    sx
  }
  paste(vapply(params, fmt_one, character(1)), collapse = "_")
}

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

save_cov = FALSE

ii = 1
igrid = c(1:3)
for(ii in igrid){
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 1
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = make_param_tag(params)
    
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
    trim_params = make_param_tag(params)
    
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
Tables[[1]]
Tables[[2]]
Tables[[3]]

dir.create("tables", showWarnings = FALSE)
for (i in seq_along(Tables)) {
  exp_name = tolower(names(experiments)[igrid[i]])
  write_latex_cov_table(
    tab = Tables[[i]],
    row_values = Mgrid[idx_rows],
    row_name = "M",
    caption = "",
    label = paste0("tab:Exp1_features_", exp_name),
    file = paste0("tables/SS_features_nfix_", exp_name, ".tex")
  )
}

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
    trim_params = make_param_tag(params)
    
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
  jj = 1
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = make_param_tag(params)
    
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
    trim_params = make_param_tag(params)
    
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
    trim_params = make_param_tag(params)
    
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

dir.create("tables", showWarnings = FALSE)
for (i in seq_along(Tables)) {
  exp_name = tolower(names(experiments)[i])
  write_latex_cov_table(
    tab = Tables[[i]],
    row_values = Ngrid[idx_rows],
    row_name = "n",
    caption = "",
    label = paste0("tab:Exp1_features_Mfix_", exp_name),
    file = paste0("tables/SS_features_Mfix_", exp_name, ".tex")
  )
}

## CI length ---------------------------------------------------------------
ii = 1
for(ii in 1:length(experiments)){
  
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 2
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = make_param_tag(params)
    
    # Load
    filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
    load(filename)
    
    # Plot
    oracle = sapply(ExpRes_list, function(x) quantile(x[,1], 1-alpha))
    ExpRes_qnt = lapply(ExpRes_list, function(x) apply(x[,c(2:6)],2,quantile,probs = c(0.025,0.5,0.975)) )
    ExpRes_qnt <- simplify2array(ExpRes_qnt) # 3 x 6 x LNgrid
    
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
  jj = 2
  for(jj in 1:Ncases){
    cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
    params = experiments[[ii]][[jj]]
    trim_params = make_param_tag(params)
    
    # Load
    filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
    load(filename)
    
    # Plot
    ParEst_mcmc = lapply(ExpRes_list, function(x) x[,c(7:16)] ) # Kn in included here
    ParEst_mcmc <- simplify2array(ParEst_mcmc) # Nrep x Nparams x LNgrid
    params_names = c(colnames(ExpRes_list[[1]])[7:16])
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








## Brutta ------------------------------------------------------------------
n = 2000
ii = 1
jj = 1
idx_keep = 3:5
Mgrid_brutta = seq(100, 10000, by = 500)
Mvals_brutta = Mgrid_brutta[idx_keep]
save_name_base = paste0("save/SS_features_nfix_")

name = names(experiments)[ii]
Ncases = length(experiments[[ii]])
cat("\n ---- ",name," ",jj,"/",Ncases," ---- \n")
params = experiments[[ii]][[jj]]
trim_params = make_param_tag(params)

# Load
filename = paste0(save_name_base,name,"_",trim_params,".Rdat")
load(filename)

Brutta_FB = lapply(idx_keep, function(idx) {
  x = ExpRes_list[[idx]]
  M_val = Mgrid_brutta[idx]
  data.frame(
    M = M_val,
    Kn = x[,16],
    Mstar = as.integer(M_val-x[,16]),
    p_Mstar_not0 = as.integer(x[,16]<M_val),
    UB_FB = x[,6],
    a_FB = x[,14],
    b_FB = x[,15],
    Mmax = x[,1]
  )
})
Brutta_MBP = lapply(idx_keep, function(idx) {
  x = ExpRes_list[[idx]]
  a = x[,10]; b = x[,11]; m = x[,12]; q = x[,13]; Kn = x[,16]
  kappa_n = exp( lgamma(b+n)-lgamma(b)+lgamma(a+b)-lgamma(a+b+n) )
  
  data.frame(
    M = Mgrid_brutta[idx],
    Kn = x[,16],
    Mstar = (m+Kn) * q*kappa_n/(1 - 1*kappa_n),
    UB_MBP = x[,5],
    a_MBP = x[,10],
    b_MBP = x[,11],
    m_MBP = x[,12],
    q_MBP = x[,13],
    Mmax = x[,1]
  )
})


Tables[[1]][3:5,10:15]
lapply(Brutta_FB, function(x) format(round(colMeans(x), 4), scientific = FALSE))
lapply(Brutta_MBP, function(x) format(round(colMeans(x), 4), scientific = FALSE))



M_Kn_Mstar_1 = Brutta_FB[[1]][,c(1:3)]

kappa_est = apply(Brutta_FB[[1]][1:10,c(6,7)],1,sum)+1
Betamu_est = apply(Brutta_FB[[1]][1:10,c(6,7)],1,function(y) y[1]/sum(y))


Rmax = 100
n = 2000
M = 1100
Kn_grid_plot = c(1097,1098,1099)
alpha <- alfa <- 0.05
tol = sqrt(.Machine$double.eps) 

mean_grid = seq(1e-5,1-1e-5,length.out = 1000)
kappa_grid = c(kappa_est,100,150,200) 
mycol_ub = c(
  gray.colors(length(kappa_est), start = 0.2, end = 0.7),
  c("#96D84B", "#CDE030", "#FDE333")
)
xmax = max(mean_grid); xmin = min(mean_grid)
xpos = seq(xmin,xmax,length.out = 10)
xlabs = as.character(floor(xpos*100))
xlim_plot = c(xmin,xmax)

par( mfrow = c(1,3), mar = c(3.5,6,1.5,1), mgp=c(4.5,1,0), bty = "l", las = 1, cex.lab = 2 )
for(Kn in Kn_grid_plot){
  Mmax_selected = Brutta_FB[[1]][Brutta_FB[[1]][,2] == Kn,8]
  UB_mat = matrix(-1, nrow = length(mean_grid), ncol = length(kappa_grid))
  
  for(i in seq_along(mean_grid)){
    vmax = mean_grid[i] * (1 - mean_grid[i])
    for(j in seq_along(kappa_grid)){
      v = vmax/kappa_grid[j] 
      ab = compute_ab_beta(mean_grid[i],v)
      UB_mat[i,j] = min(1, exp(compute_log_UBMarkov_FB(Rmax,ab[1],ab[2],n,Kn,M,alpha) ) )
    }
  }
  
  ymax = (11/10) * max(UB_mat); ymin = (10/11) * min(UB_mat)
  ypos = seq(ymin,ymax,length.out = 5)
  ylabs = as.character(round(ypos,2))
  
  plot(0,0,  yaxt = "n", xaxt = "n",
       xlab = "", ylab = "bound",
       xlim = c(0,0.1) , ylim = c(0,0.01), 
       main = paste0("Kn = ", Kn),
       type = "n")
  grid(lty = 1,lwd = 1, col = "gray90" )
  axis(side = 2, at = ypos, labels = ylabs, cex.axis = cex.axis )
  axis(side = 1, at = xpos, labels = xlabs, cex.axis = cex.axis )
  mtext("mean * 100", side = 1, line = 2.5, cex = cex.axis)
  for(j in seq_along(kappa_grid)){
    points( x = mean_grid, y = UB_mat[,j], 
            type = "l", lwd = 5, col = mycol_ub[j] )
  }
  abline(h = Mmax_selected, lty = 2, col = "red", lwd = 0.5)
  abline(v = Betamu_est, lty = 2, col = "black", lwd = 0.5)
  # legend("topright",
  #        legend = paste0("kappa = ", format(kappa_grid, scientific = FALSE)),
  #        col = mycol_ub, lwd = 5,
  #        bty = "n", cex = cex.legend)
}
