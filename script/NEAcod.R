library(stockassessment)
source("script/mohnsDist.R")

#Read data
cn<-read.ices("data/NEACod/cn.dat")
cw<-read.ices("data/NEACod/cw.dat")
dw<-read.ices("data/NEACod/dw.dat")
lw<-read.ices("data/NEACod/lw.dat")
mo<-read.ices("data/NEACod/mo.dat")
nm<-read.ices("data/NEACod/nm.dat")
pf<-read.ices("data/NEACod/pf.dat")
pm<-read.ices("data/NEACod/pm.dat")
sw<-read.ices("data/NEACod/sw.dat")
lf<-read.ices("data/NEACod/lf.dat")
surveys<-read.ices("data/NEACod/survey.dat")
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
conf = loadConf(dat,file = "data/NEACod/conf.cfg")
par = defpar(dat,conf)
fit = sam.fit(dat,conf,par)
nya = 5
ret = retro(fit,year = nya,ncores = nya)

nsim = 1020
nsimSave = 1000
ncores = 5

#Distribution in official assessment
start_time <- Sys.time()
rets = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years )
end_time <- Sys.time()
end_time - start_time

#Plot (1-p)*100 % confidence intervals
plotMohnDist(rets,ret, examples = FALSE, p = 0.05)
#Plot realizations in simulation experiment
plotMohnDist(rets,ret, examples = TRUE, xlim = c(2010,2020))

#Save sample of Mohn's rho
saveMohnsSamples(rets,file = "results/NEACodOfficial.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(rets)/nsim))

#Retrospective plot
png(file = "figures/retroSSBNEACod.png")
rr= rets[1:nsimSave]
qq = quantileRetro(rr)
lQ = qq$lQ
uQ = qq$uQ
ssbplot(ret,xlim =c(2000,2020), main = "SSB retro plot NEA cod")
years = fit$data$years
library(plotrix)
for(y in 1:nya){
  rl = as.list(fit$sdrep,what = "Est",report = TRUE)
  est = exp(rl$logssb[length(years)-y])

  plotCI(max(years)-y, est, ui = est + est*uQ$SSB[y],li = est+ est*lQ$SSB[y], add = TRUE)
}
dev.off()


#Scenarios with assumtion violations. Run the following lines, and then the script "potsStrength.R"
nsim = 210
nsimSave = 200
#M*3
scaleM = seq(1,3,length.out = nya+1)
retsM3 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                   scaleM = scaleM )
saveMohnsSamples(retsM3,file = "results/NEACodM3.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsM3)/nsim))


#M/3
scaleM = seq(1,1/3,length.out = nya+1)
retsM05 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                    scaleM = scaleM)
saveMohnsSamples(retsM05,file = "results/NEACodM05.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsM05)/nsim))

#Q*3
scaleQ = seq(1,3,length.out = nya+1)
retsQ3 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                   scaleQ = scaleQ)
saveMohnsSamples(retsQ3,file = "results/NEACodQ3.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsQ3)/nsim))


#C/3
scaleC = seq(1,1/3,length.out = nya+1)
retsC05 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                    scaleC = scaleC)
saveMohnsSamples(retsC05,file = "results/NEACodC05.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsC05)/nsim))


#varF*3
scaleSdF = seq(1,sqrt(3),length.out = nya+1)
retsSdF3 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                      scaleSdF = scaleSdF )
saveMohnsSamples(retsSdF3,file = "results/NEACodSdF3.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsSdF3)/nsim))


#var obsC *3
scaleSdObsC= seq(1,sqrt(3),length.out = nya+1)
retsSdObsC3 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                        scaleSdObsC =  scaleSdObsC)
saveMohnsSamples(retsSdObsC3,file = "results/NEACodSdObsC3.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsSdObsC3)/nsim))


#var obsI *3
scaleSdObsI= seq(1,sqrt(3),length.out = nya+1)
retsSdObsI3 = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                        scaleSdObsI =  scaleSdObsI)
saveMohnsSamples(retsSdObsI3,file = "results/NEACodSdObsI3.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsSdObsI3)/nsim))


#Remove corObs
scaleCorObs= TRUE
retsCorObs = mohnsDist(ret,seed = 123,nsim = nsim,ncores = ncores, rec.years = ret[[nya]]$data$years,
                       scaleCorObs =  scaleCorObs)
saveMohnsSamples(retsCorObs,file = "results/NEACodCorObs.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsCorObs)/nsim))



############################With AR in F
library(stockassessment)
source("script/mohnsDist.R")
cn<-read.ices("data/NEACod/cn.dat")
cw<-read.ices("data/NEACod/cw.dat")
dw<-read.ices("data/NEACod/dw.dat")
lw<-read.ices("data/NEACod/lw.dat")
mo<-read.ices("data/NEACod/mo.dat")
nm<-read.ices("data/NEACod/nm.dat")
pf<-read.ices("data/NEACod/pf.dat")
pm<-read.ices("data/NEACod/pm.dat")
sw<-read.ices("data/NEACod/sw.dat")
lf<-read.ices("data/NEACod/lf.dat")
surveys<-read.ices("data/NEACod/survey.dat")
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
conf = loadConf(dat,file = "data/NEACod/conf.cfg")
conf$corFlag=2
par = defpar(dat,conf)
fitAR1 = sam.fit(dat,conf,par)
nya = 5
retAR1 = retro(fitAR1,year = nya,ncores = nya)

#Simulate distribution of Mohn's rho
nsim = 1020
nsimSave = 1000
ncores = 10
retsAR1 = mohnsDist(retAR1,seed = 123,nsim = nsim,ncores = ncores, rec.years = retAR1[[nya]]$data$years )

#Plot distribution of Mohn's rho
plotMohnDist(retsAR1,retAR1, examples = FALSE, p = 0.05)
#Plot realizations in simulation experiment
plotMohnDist(retsAR1,retAR1, examples = TRUE, xlim = c(2010,2020), p = 0.05)

#Save Mohn's rho samples
saveMohnsSamples(retsAR1,file = "results/NEACodAR1.txt",nsim = nsimSave)
print(paste0("proportion converged: ", length(retsAR1)/nsim))

#Retrospective plots
png(file = "figures/retroFNEACod.png")
rets= retsAR1[1:nsimSave]
qq = quantileRetro(rets)
lQ = qq$lQ
uQ = qq$uQ
fbarplot(retAR1,xlim =c(2000,2020), main = "F-bar retro plot NEA cod",cex.lab= 1.1,cex.main = 1.5)
years = fitAR1$data$years-1
library(plotrix)
for(y in 1:nya){
  rl =as.list(fitAR1$sdrep,what = "Est", report = TRUE)
  est = exp(rl$logfbar)[length(years)-y-1]
  plotCI(max(years)-y, est, ui = est + est*uQ$F[y],li = est+ est*lQ$F[y], add = TRUE)
}
dev.off()


