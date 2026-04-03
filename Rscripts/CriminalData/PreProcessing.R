# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"bnp_upperbounds/Rscripts/CriminalData/")
setwd(wd)
# Load data ---------------------------------------------------------------
A <- read.csv(file="NDRANGHETAMAFIA_2M.csv",header=TRUE, stringsAsFactors = TRUE)
# A is 156x48
n_people = nrow(A) # 156
n_meetings = ncol(A)-1 # 47, first column reports the names

# Colors ------------------------------------------------------------------
save_img = FALSE
width = 12; height = 6
bw_col = c("grey90","black")


# Functions ---------------------------------------------------------------
source("../../R/Rfunctions.R")

# Pre processing ----------------------------------------------------------
# missed meetings
A[23,20] <- 1; A[48,23] <- 1; A[78,35] <- 1; A[115,22] <- 1

# suspects who never attended a summit
sel_empty <- which(apply(as.matrix(A[,-1]),1,sum)==0)

# suspects who have not been recognized during the investigation process (we only have the name, not the surname)
A[c(38,105,106,125,135),1]

# indicators of the two groups of suspects to be excluded
sel_empty <- c(c(sel_empty),c(38,105,106,125,135))

# remove these suspects from the dataset
A <- A[-sel_empty,]

# create a vector with the actors' names
actors <- A[,1]
actors <- droplevels(actors)


## Plot raw data -----------------------------------------------------------

A <- as.matrix(A[,-1])
n_people = nrow(A) # 146
n_meetings = ncol(A) # 47

par(mfrow = c(1,1), mar = c(1,3,3,1), mgp=c(1.5,0.5,0))
image( 1:n_people, 1:n_meetings,
       A[,ncol(A):1],
       zlim = c(0,1),
       col = bw_col,
       ylab = "Meetings",
       xlab = " ",
       main = " ",
       axes = FALSE ); box()
axis(2, at = 1:n_meetings,
     labels = n_meetings:1,
     cex.axis = 0.7 )
axis(3, at = 1:n_people, 
     labels = 1:n_people,
     cex.axis = 0.7); mtext("People", side = 3, line = 2)


## Plot network data -----------------------------------------------------------

Anet <- A%*%t(A) # (146x146 matrix)
Anet <- (Anet>0)*1
diag(Anet) <- 0
rownames(Anet) <- colnames(Anet) <- c(1:dim(Anet)[1])

par(mfrow = c(1,1), mar = c(1,3,3,1), mgp=c(1.5,0.5,0))
image( 1:n_people, 1:n_people,
       Anet[,ncol(Anet):1],
       zlim = c(0,1),
       col = bw_col,
       ylab = "People",
       xlab = " ",
       main = " ",
       axes = FALSE ); box()
axis(2, at = 1:n_people,
     labels = n_people:1,
     cex.axis = 0.7 )
axis(3, at = 1:n_people, 
     labels = 1:n_people,
     cex.axis = 0.7); mtext("People", side = 3, line = 2)



# Add Locale/Role information --------------------------------------------------

# Locale membership ("OUT": Suspects not belonging to La Lombardia. 
#                   "MISS": Information not available)

Locale <- c("C","OUT","A","MISS","O","A","MISS","D","D","D","D","D","C","P","L","L","Q","MISS","B","OUT","B","B","I","MISS","OUT","D","A","O","N","N","H","OUT","D","E","G","G","L","A","OUT","Q","C","OUT","Q","L","C","MISS","C","C","F","C","OUT","D","A","B","B","E","M","MISS","C","C","C","B","H","C","C","E","E","E","E","C","MISS","L","A","A","E","E","C","E","E","E","C","MISS","OUT","C","C","E","G","A","A","B","I","I","A","B","B","OUT","I","A","G","N","E","D","F","OUT","OUT","C","D","C","MISS","MISS","C","MISS","E","E","C","MISS","OUT","B","L","A","D","D","O","MISS","B","D","O","D","D","A","A","I","C","MISS","MISS","MISS","A","A","F","E","C","Q","H","B","B","B")   
length(which(Locale == "MISS"))
# Leadership role ("miss": Information not available)

Role <- c("aff","aff","aff","miss","aff","boss","miss","boss","boss","aff","aff","aff","aff","aff","aff","aff","aff","miss","aff","boss","boss","boss","boss","miss","boss","boss","aff","aff","aff","boss","aff","boss","aff","aff","aff","aff","aff","aff","boss","aff","aff","aff","boss","aff","aff","miss","aff","aff","aff","aff","boss","aff","aff","aff","aff","aff","boss","miss","aff","aff","aff","aff","boss","boss","boss","aff","boss","aff","aff","boss","miss","aff","aff","boss","boss","aff","aff","aff","aff","aff","aff","miss","boss","aff","aff","aff","aff","aff","aff","boss","aff","boss","aff","aff","aff","aff","boss","boss","boss","boss","aff","aff","aff","aff","aff","aff","aff","boss","miss","miss","aff","miss","aff","aff","aff","miss","aff","aff","boss","aff","aff","aff","aff","miss","aff","aff","boss","aff","boss","aff","aff","aff","aff","miss","miss","miss","aff","aff","boss","aff","aff","aff","boss","aff","aff","boss")

# indicators of those suspects who are not known to be part of the organization La Lombardia 
sel_miss <- which(Locale=="MISS" | Locale=="OUT")

# remove such suspects from the dataset
Locale_temp <- Locale[-sel_miss]
Role_temp <- Role[-sel_miss]
actors_temp <- actors[-sel_miss]

A <- A[-sel_miss,] # (118x47)
Anet <- Anet[-sel_miss,-sel_miss] #(118x118)
rownames(Anet) <- colnames(Anet) <- c(1:dim(Anet)[1]) 

# interaction between Locale and Role
RoleLocale_temp <- paste(Role_temp,Locale_temp,sep="_")


## Plot Locali abundances --------------------------------------------------

Nj_locale = table(Locale_temp)
Nj_locale = sort(Nj_locale, decreasing = TRUE)
bp1 = barplot(Nj_locale)


if(save_img)
  pdf("img/Barplot_AbundanceLocale.pdf", width=width, height=height)
par(mfrow = c(1,1), mgp=c(1.5,0.5,0), mar = c(2,2.5,1,0),cex = 2)
barplot( height = Nj_locale, 
         names.arg = "", las = 2, col = "darkred", border = NA,
         main = " ", ylab = "#obs.", yaxt = "n" )
axis( side = 2, at = seq(0, max(Nj_locale), by=5), las = 1)
text( x = bp1, y = -1, 
      labels = names(Nj_locale), cex = 0.9,
      srt = 0, adj = 0.5, xpd = TRUE)
if(save_img)
  dev.off()


## Plot raw data - reduced -----------------------------------------------------------

n_people = nrow(A) # 118
n_meetings = ncol(A) # 47

par(mfrow = c(1,1), mar = c(1,3,3,1), mgp=c(1.5,0.5,0))
image( 1:n_people, 1:n_meetings,
       A[,ncol(A):1],
       zlim = c(0,1),
       col = bw_col,
       ylab = "Meetings",
       xlab = " ",
       main = " ",
       axes = FALSE ); box()
axis(2, at = 1:n_meetings,
     labels = n_meetings:1,
     cex.axis = 0.7 )
axis(3, at = 1:n_people, 
     labels = 1:n_people,
     cex.axis = 0.7); mtext("People", side = 3, line = 2)


## Plot network data -----------------------------------------------------------

par(mfrow = c(1,1), mar = c(1,3,3,1), mgp=c(1.5,0.5,0))
image( 1:n_people, 1:n_people,
       Anet[,ncol(Anet):1],
       zlim = c(0,1),
       col = bw_col,
       ylab = "People",
       xlab = " ",
       main = " ",
       axes = FALSE ); box()
axis(2, at = 1:n_people,
     labels = n_people:1,
     cex.axis = 0.7 )
axis(3, at = 1:n_people, 
     labels = 1:n_people,
     cex.axis = 0.7); mtext("People", side = 3, line = 2)

## Plot Binomial counts -----------------------------------------------------------
### Mod2 - Inc.4.meetings -----------------------------------------------------------
n = nrow(A)
Kn = ncol(A)
N_j = colSums(A)
N_j = sort(N_j, decreasing = TRUE)
names(N_j) = sapply(N_j, function(x) paste0("M",x))
bp1 = barplot(N_j)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2.5,3.5,1,0))
barplot( height = N_j, 
         names.arg = "", las = 2, col = "darkred", border = NA,
         main = " ", ylab = "#obs.", yaxt = "n" )
axis( side = 2, at = seq(0, max(N_j), by=2), las = 1)
text( x = bp1, y = -1, 
      labels = names(N_j), cex = 0.2,
      srt = 0, adj = 0.5, xpd = TRUE)


### Mod2 - Inc.4.people -----------------------------------------------------------
n = ncol(A)
Kn = nrow(A)
N_j = rowSums(A)
N_j = sort(N_j, decreasing = TRUE)
names(N_j) = sapply(N_j, function(x) paste0(""))
bp1 = barplot(N_j)

par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2.5,3.5,1,0))
barplot( height = N_j, 
         names.arg = "", las = 2, col = "darkred", border = NA,
         main = " ", ylab = "#obs.", yaxt = "n" )
axis( side = 2, at = seq(0, max(N_j), by=2), las = 1)
text( x = bp1, y = -1, 
      labels = names(N_j), cex = 0.2,
      srt = 0, adj = 0.5, xpd = TRUE)


ord <- order(N_j, decreasing = TRUE)
N_j_sorted <- N_j[ord]
Role_sorted <- Role_temp[ord]

cols <- ifelse(Role_sorted == "boss", "darkred", "steelblue") # colors by role

if(save_img)
  pdf("img/Mod3_Inc4people_barplot.pdf", width = width, height = height)
par(mfrow = c(1,1), mgp=c(1,0.5,0), mar = c(1,2,1,0), cex = 2)
bp1 <- barplot(height = N_j_sorted,
               names.arg = "",
               col = cols, border = NA,
               main = " ", ylab = "", yaxt = "n")
axis(2, at = seq(0, max(N_j_sorted), by = 2), las = 1)
legend("topright", legend = c("boss", "affil"),
       fill = c("darkred", "steelblue"), border = NA, bty = "n")
if(save_img)
  dev.off()
# Accumulation curves -----------------------------------------------------
seed = 42
set.seed(seed)


## Local abundances ------------------------------------------------------
n_reorderings = 50
orderings = lapply(1:n_reorderings, function(j) sample(1:n_people, size = n_people))
AccCurves_list = lapply(orderings, function(o) Locale_temp[o])
AccCurves_list = lapply(AccCurves_list, function(x){
  sapply(1:length(x), function(j) { y = x[1:j]; length(table(y)) })
})
AccCurves_mat <- do.call(rbind, AccCurves_list) #(n_reorderings x n_people)
AccCurves_res = apply(AccCurves_mat, 2, quantile, probs = c(0.025,0.5,0.975))

# Plot
ngrid = 1:n_people 

if(save_img)
  pdf("img/AccCurves_LocaleAbundances.pdf", width=width, height=height)
par(mfrow = c(1,1), mar = c(2.5,2.5,0,1), mgp=c(1.5,0.5,0), bty = "l", cex = 2)
plot(x = 0, y = 0, type = "n",
     main = " ", xlab = "n", ylab = "#Locali",
     ylim = c(0,max(AccCurves_res)+1),
     xlim = c(0,n_people+1),
     pch = 1) # init plot
polygon( c(ngrid, rev(ngrid)),
         c(AccCurves_res[1,], rev(AccCurves_res[3,])),
         col = "grey75",
         border = NA) # plot in-sample bands
points(x = ngrid, y = AccCurves_res[2,], type = "l", lwd = 3) # plot mean obs
if(save_img)
  dev.off()

## Meetings level ------------------------------------------------------
seed = 12131
n_reorderings = 50
AccCurves_mod2 = rarefaction.array(object = A, n_reorderings = n_reorderings, seed = seed)


ngrid = 1:n_people 
par(mfrow = c(1,1), mar = c(3,3,2,1), mgp=c(1.5,0.5,0), bty = "l")
plot(x = 0, y = 0, type = "n",
     main = " ", xlab = "#obs.", ylab = "#Meetings",
     ylim = c(0,max(AccCurves_mod2)+1),
     xlim = c(0,n_people+1),
     pch = 1) # init plot
polygon( c(ngrid, rev(ngrid)),
         c(AccCurves_mod2[1,], rev(AccCurves_mod2[3,])),
         col = "grey75",
         border = NA) # plot in-sample bands
points(x = ngrid, y = AccCurves_mod2[2,], type = "l", lwd = 3) # plot mean obs


## People level ------------------------------------------------------
seed = 12131
n_reorderings = 50
AccCurves_mod2 = rarefaction.array(object = t(A), n_reorderings = n_reorderings, seed = seed)


ngrid = 1:n_meetings 
par(mfrow = c(1,1), mar = c(3,3,2,1), mgp=c(1.5,0.5,0), bty = "l")
plot(x = 0, y = 0, type = "n",
     main = " ", xlab = "#obs.", ylab = "#People",
     ylim = c(0,max(AccCurves_mod2)+1),
     xlim = c(0,n_meetings+1),
     pch = 1) # init plot
polygon( c(ngrid, rev(ngrid)),
         c(AccCurves_mod2[1,], rev(AccCurves_mod2[3,])),
         col = "grey75",
         border = NA) # plot in-sample bands
points(x = ngrid, y = AccCurves_mod2[2,], type = "l", lwd = 3) # plot mean obs










