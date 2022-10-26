#devtools::install_github("fishfollower/SAM/stockassessment",ref = "devMinor") #Minor correction of how time in year comes into the biomass index
library(stockassessment)
source("script/mohnsDist.R")

#Read data
cn<-read.ices("data/coastalCod/cn.dat")
cw<-read.ices("data/coastalCod/cw.dat")
dw<-read.ices("data/coastalCod/dw.dat")
lw<-read.ices("data/coastalCod/lw.dat")
mo<-read.ices("data/coastalCod/mo.dat")
nm<-read.ices("data/coastalCod/nm.dat")
pf<-read.ices("data/coastalCod/pf.dat")
pm<-read.ices("data/coastalCod/pm.dat")
sw<-read.ices("data/coastalCod/sw.dat")
lf<-read.ices("data/coastalCod/lf.dat")
surveys<-read.ices("data/coastalCod/survey.dat")
dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn,
                    prop.mature=mo,
                    stock.mean.weight=sw,
                    catch.mean.weight=cw,
                    dis.mean.weight=dw,
                    land.mean.weight=lw,
                    prop.f=pf,
                    prop.m=pm,
                    natural.mortality=nm,
                    land.frac=lf)



#Fit model and calculate retro
conf = loadConf(dat,file = "data/coastalCod/conf.cfg")
par = defpar(dat,conf)
fit = sam.fit(dat,conf,par)
nya = 5
ret = retro(fit,year = nya,ncores = nya)

#Calculate distribution of Mohn's rho
nsim = 1100
nsimSave = 1000 #Simulation to save
ncores = 5

start_time <- Sys.time()
rets = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years )
end_time <- Sys.time()
end_time - start_time

#Plot (1-p)*100 % confidence intervals of Mohn's rho
plotMohnDist(rets,ret, examples = FALSE, p = 0.05)

#Plot reaslisations in simulation experiment
plotMohnDist(rets,ret, examples = TRUE, xlim = c(2010,2020))

#Save simulated Mohn's rho
saveMohnsSamples(rets,file = "results/costalCodOfficial.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(rets)/nsim))

#Retrospective plots
png(file = "figures/retroSSBNCCCod.png")
rr= rets[1:nsimSave]
qq = quantileRetro(rr)
lQ = qq$lQ
uQ = qq$uQ
ssbplot(ret,xlim =c(2000,2020), main = "SSB retro plot NCC cod")
years = fit$data$years
library(plotrix)
for(y in 1:nya){
    rl = as.list(fit$sdrep,what = "Est",report = TRUE)
    est = exp(rl$logssb[length(years)-y])

    plotCI(max(years)-y, est, ui = est + est*uQ$SSB[y],li = est+ est*lQ$SSB[y], add = TRUE)
}
dev.off()

png(file = "figures/retroFNCCCod.png")
rr= rets[1:nsimSave]
qq = quantileRetro(rr)
lQ = qq$lQ
uQ = qq$uQ
fbarplot(ret,xlim =c(2000,2020), main = "F-bar retro plot NCC cod")
years = fit$data$years
library(plotrix)
for(y in 1:nya){
  rl = as.list(fit$sdrep,what = "Est",report = TRUE)
  est = exp(rl$logfbar[length(years)-y])

  plotCI(max(years)-y, est, ui = est + est*uQ$F[y],li = est+ est*lQ$F[y], add = TRUE)
}
dev.off()


#Scenarios with assumption violations. Run the following lines to save simulations, and then the script "potsStrength.R"
nsim = 250
nsimSave = 200

#M3
scaleM = seq(1,3,length.out = nya+1)
retsM3 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                   scaleM = scaleM)
saveMohnsSamples(retsM3,file = "results/costalCodM3.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsM3)/nsim))


#M/3
scaleM = seq(1,1/3,length.out = nya+1)
retsM05 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                    scaleM = scaleM)
saveMohnsSamples(retsM05,file = "results/costalCodM05.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsM05)/nsim))

#Q3
scaleQ = seq(1,3,length.out = nya+1)
retsQ3 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                   scaleQ = scaleQ )
saveMohnsSamples(retsQ3,file = "results/costalCodQ3.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsQ3)/nsim))

#C/3
scaleC = seq(1,1/3,length.out = nya+1)
retsC05 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                    scaleC = scaleC )
saveMohnsSamples(retsC05,file = "results/costalCodC05.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsC05)/nsim))

#varF *3
scaleSdF = seq(1,sqrt(3),length.out = nya+1)
retsSdF3 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                      scaleSdF = scaleSdF )
saveMohnsSamples(retsSdF3,file = "results/costalCodSdF3.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsSdF3)/nsim))


#Var obsC *3
scaleSdObsC= seq(1,sqrt(3),length.out = nya+1)
retsSdObsC3 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                       scaleSdObsC =  scaleSdObsC)
saveMohnsSamples(retsSdObsC3,file = "results/costalCodSdObsC3.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsSdObsC3)/nsim))


#Var obsI *3
scaleSdObsI= seq(1,sqrt(3),length.out = nya+1)
retsSdObsI3 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                        scaleSdObsI =  scaleSdObsI)
saveMohnsSamples(retsSdObsI3,file = "results/costalCodSdObsI3.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsSdObsI3)/nsim))


#Remove corObs
scaleCorObs= TRUE
retsCorObs = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                        scaleCorObs =  scaleCorObs)
saveMohnsSamples(retsCorObs,file = "results/costalCodCorObs.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsCorObs)/nsim))




