

# Brutta ------------------------------------------------------------------
n = 1000
Nrep = 5
Mmax = 10000
sigma = 0.25; theta = 10
gamma = 0.15; Lambda = 100; M = 102
AccCrv_length = 20
seed = 23131


params = c(sigma,theta); model = "PD"
ptrue1 = SimModel_generic(model,Mmax,Nrep,params)
colSums(ptrue1)

gof_pd = GOF_generic(model,n,Mmax,Nrep,params,seed,AccCrv_length)
gof_pd$Envelop_qnt[,1:10]
gof_pd$Freq.Rare
gof_pd$AccCrv[,1:10]

params = c(gamma,Lambda); model = "FDP"
ptrue2 = SimModel_generic(model,Mmax,Nrep,params)
colSums(ptrue2)

gof_FDP = GOF_generic(model,n,Mmax,Nrep,params,seed,AccCrv_length)
gof_FDP$Envelop_qnt[,1:10]
gof_FDP$Freq.Rare
gof_FDP$AccCrv[,1:10]


params = c(gamma,M); model = "DirMulti"
ptrue3 = SimModel_generic(model,Mmax,Nrep,params)
colSums(ptrue3)

gof_DM = GOF_generic(model,n,Mmax,Nrep,params,seed,AccCrv_length)
gof_DM$Envelop_qnt[,1:10]
gof_DM$Freq.Rare
gof_DM$AccCrv[,1:10]

ptrue_mat = ptrue3
