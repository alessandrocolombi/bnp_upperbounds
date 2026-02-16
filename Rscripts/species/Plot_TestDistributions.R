params_zipfs = list(0.9,1.02,2)
params_geom = list(0.85,0.9,0.95)
params_unif = list(NA)
params_negbin = list(c(1,0.003),c(5,0.003),c(1,0.01))
experiments = list("Zipfs" = params_zipfs,
                   "Geom" = params_geom,
                   "Uniform" = params_unif,
                   "NegBin" = params_negbin)


M = 1000
Mplot = 100
n = 500

ii = 2
for(ii in 1:length(experiments)){
  name = names(experiments)[ii]
  Ncases = length(experiments[[ii]])
  jj = 1
  for(jj in 1:Ncases){
    Nparams = length(experiments[[ii]][[jj]])
    params = experiments[[ii]][[jj]]
    ptrue = sim_generic_species(name,M,params)
    ptrue = sort(ptrue, decreasing = TRUE)
    par( mfrow = c(1,1), mar = c(3.5,4.25,1,1), mgp=c(2.75,1,0), bty = "l", las = 1, cex.lab = 1 )
    plot(x = 1:Mplot, y = ptrue[1:Mplot], type = "p", pch = 16, cex = 0.5,
         main = paste0(name,", ",ii,", ",jj), xlab = "", ylab = "")
  }
}

