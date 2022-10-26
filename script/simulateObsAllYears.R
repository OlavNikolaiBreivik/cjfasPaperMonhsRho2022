library(stockassessment)#Obs, need tiny modification of SAM such that it adreport N and F in first year which is used in the forecast
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
fitOriginal = sam.fit(dat,conf,par)
retOriginal = retro(fitOriginal,year = 5,ncores = 5)

simFullData = function(fit,nsim = 50,seed = NULL, processNoiseF = TRUE,ncores = 1,rec.years = NA){
  nya = length(fit$data$years)
  if(!is.null(seed)){
    set.seed(seed)
  }

  simDataList<-list()
  simN = list()
  simF = list()

  nm <- fit$data$natMor
  cw <- fit$data$catchMeanWeight
  dw <- fit$data$disMeanWeight
  lw <- fit$data$landMeanWeight
  lf <- fit$data$landFrac
  mo <- fit$data$propMat
  sw <- fit$data$stockMeanWeight
  pf <- fit$data$propF
  pm <- fit$data$propM
  nn = dim(nm)[2]

  scaleM = rep(1,99)
  scaleSdF = rep(1,99)
  scaleSdObsC = rep(1,99)
  scaleC = rep(1,99)
  scaleQ = rep(1,99)
  whichQScale = 1:99
  scaleSdObsI = rep(1,99)
  scaleSdObsC = rep(1,99)
  corF = FALSE
  scaleCorObs = FALSE

  if(is.na(rec.years[[1]])){
    rec.years  =fit$data$years
  }

  processNoiseF = TRUE
  fc <- forecast(fit, fscale=rep(1,nya), nosim=nsim, processNoiseF=processNoiseF, savesim = TRUE,rec.years = rec.years,
                 scaleM = scaleM,scaleSdF = scaleSdF,corF = corF, year.base = min(fit$data$years))

  ssbplot(fc)

  for(j in 1:nsim){
    cnnew <- getFleet(fit,1)[1,,drop = FALSE]
    surveysnew <- list()
    for(ff in 2:fit$data$noFleets){
      surveysnew[[ff-1]] <- getFleet(fit,ff)[1,,drop = FALSE]
      attr(surveysnew[[ff-1]],"time") <- rep(fit$data$sampleTimes[ff],2)
    }

    simN[[j]] = matrix(0,nrow = nn, ncol = nya)
    simF[[j]] = matrix(0,nrow = nn, ncol = nya)
    rownames(simF[[j]]) = fit$conf$minAge:fit$conf$maxAge
    rownames(simN[[j]]) = fit$conf$minAge:fit$conf$maxAge
    simN[[j]][,1] = exp(fc[[1]]$sim[j,1:nn])
    simF[[j]][,1] = exp(fc[[1]]$sim[j,attr(fc,"fit")$conf$keyLogFsta[1,]+(nn+1)])
    for(y in 1:nya){
      nmThis = nm[y,]*scaleM[y]
      swThis = sw[y,]
      moThis = mo[y,]

      theFleets<-list()
      N<-exp(fc[[y]]$sim[j,1:nn])
      F<-exp(fc[[y]]$sim[j,attr(fc,"fit")$conf$keyLogFsta[1,]+(nn+1)])

      simN[[j]][,y] = N
      simF[[j]][,y] = F
      Cpred <- F/(F+nmThis)*(1-exp(-F-nmThis))*N


      covCatch <- attr(fc,"fit")$rep$obsCov[[1]]
      Csim <- exp(log(Cpred)+rmvnorm(1,rep(0,nrow(covCatch)),covCatch*scaleSdObsC[y]^2))*scaleC[y]


      theFleets[[1]] <- Csim

      f <- Vectorize(function(ii)ifelse(ii==-1,NA,attr(fc,"fit")$opt$par[ii+1]))
      Q <- t(exp(apply(attr(fc,"fit")$conf$keyLogFpar,1,f)))

      Qtmp = Q #Applied when defining plus group.
      for(ff in 2:fit$data$noFleets){
        if(fit$conf$maxAgePlusGroup[ff]==1){
          Q[ff,which(is.na(Q[ff,]))] = Q[ff,max(which(Q[ff,]>0))]
        }
        if(ff %in% whichQScale){
            Q[ff,] <- Q[ff,]*scaleQ[y]
        }

      }

      tau <- attr(fc,"fit")$data$sampleTimes

      for(ff in 2:fit$data$noFleets){
        Qf <- Q[ff,]
        if(fit$data$fleetTypes[[ff]]==2){
          Sfpred <- Qf*N*exp(-tau[ff]*(F+nmThis))
          if(fit$conf$maxAgePlusGroup[ff]==1){
            Sfpred[max(which(Qtmp[ff,]>0))] = sum(Sfpred[max(which(Qtmp[ff,]>0)): nn])
          }

          if(sum(is.na(Qtmp[ff,]))>0){
            Sfpred <- Sfpred[-which(is.na(Qtmp[ff,]))]
          }

          if(!scaleCorObs){
            Sfcov <- attr(fc,"fit")$rep$obsCov[[ff]]
          }else{
            Sfcov <- attr(fc,"fit")$rep$obsCov[[ff]]
            sdCatch = sqrt(diag(Sfcov))
            Sfcov = diag(sdCatch)*diag(dim(Sfcov)[1])*diag(sdCatch)
          }
        }else if(fit$data$fleetTypes[[ff]]==3){
          if(fit$conf$keyBiomassTreat[[ff]]==0){
            #SSB-index
            Sfpred <- Qf[1]*sum(N*exp(-tau[ff]*(F+nmThis))*swThis*moThis)
          }else if(fit$conf$keyBiomassTreat[[ff]]==5){
            #TSB-index #NB!!! Can't see propM og propF being used (but applied for SSB predictions)
            Sfpred <- Qf[1]*sum(N*exp(-tau[ff]*(F+nmThis))*swThis)
          }
          else if(fit$conf$keyBiomassTreat[[ff]]==6){
            #TSN-index
            Sfpred <- Qf[1]*sum(N*exp(-tau[ff]*(F+nmThis)))
          }else{
            stop(paste0("Mohns dist not yet implemented for bimoassTreat ",fit$conf$keyBiomassTreat[[ff]]))
          }
          Sfcov = as.matrix(exp(fit$pl$logSdLogObs[fit$conf$keyVarObs[ff,1]+1])^2)
        }
        Sfsim <- exp(log(Sfpred)+rmvnorm(1,rep(0,nrow(Sfcov)),Sfcov*scaleSdObsI[y]^2))

        theFleets[[ff]] <- Sfsim
      }

      cnnew <- addrow(cnnew,theFleets[[1]])
      for(ff in 2:fit$data$noFleets){
        surveysnew[[ff-1]] <- addrow(surveysnew[[ff-1]],theFleets[[ff]])
      }
    }

    cnnew = cnnew[-1,,drop = FALSE]
    rownames(cnnew) = fit$data$years

    for(ff in 2:fit$data$noFleets){
      surveysnew[[ff-1]] = surveysnew[[ff-1]][-1,,drop = FALSE]
      rownames(surveysnew[[ff-1]]) = fit$data$years
    }

    for(ff in 2:fit$data$noFleets){
      ss = surveysnew[[ff-1]]
      ssOrig = getFleet(fit,ff)

      keep = which((rownames(ss) %in% rownames(ssOrig)))
      ss = ss[keep,,drop = FALSE]

      ss[is.na(ssOrig)] = NA
      surveysnew[[ff-1]] =ss
      attr(surveysnew[[ff-1]],"time") <- rep(fit$data$sampleTimes[ff],2)
    }

    attributes(surveysnew)$names = attributes(fit$data)$fleetNames[-1]
    simDataList[[j]] <-setup.sam.data( surveys = surveysnew , residual.fleet =cnnew , prop.mature =mo , stock.mean.weight =sw ,
                                       catch.mean.weight =cw , dis.mean.weight =dw , land.mean.weight =lw , prop.f=pf ,
                                       prop.m=pm, natural.mortality =nm , land.frac =lf)
    attributes(simDataList[[j]])$fleetNames = attributes(fit$data)$fleetNames
  }


  if(ncores>1){
    cl <- makeCluster(ncores) #set up nodes
    clusterExport(cl, varlist=c("fit","simDataList"), envir=environment())
    lib.ver <- dirname(path.package("stockassessment"))
    clusterExport(cl, varlist="lib.ver", envir=environment())
    clusterEvalQ(cl, {library(stockassessment, lib.loc=lib.ver)})
    fitsim <- parLapply(cl,simDataList, function(f) {
      conf <- fit$conf
      parnew <- defpar(f , conf)
      out <- tryCatch(
        {
          sam.fit(f ,conf , parnew, silent=TRUE, newtonsteps=0, ignore.parm.uncertainty=TRUE)
        },
        error=function(cond) {
          return(NA)
        }
      )
      out
    })
    stopCluster(cl)
  }else {
    conf <- fit$conf
    fitsim <- lapply(simDataList, function(f) {
      parnew <- defpar(f , conf)
      sam.fit(f ,conf , parnew, silent=TRUE, newtonsteps=0, ignore.parm.uncertainty=TRUE)
    })
  }


  issueConv = rep(FALSE,nsim)
  for(i in 1:length(fitsim)){
    if(!is.na(fitsim[[i]])[1]){
      if(fitsim[[i]]$opt$convergence != 0){
        issueConv[i] = TRUE
      }
    }else{
      issueConv[i] = TRUE
    }
  }

  fitsim = fitsim[which(!issueConv)]
  simN = simN[which(!issueConv)]
  simF = simF[which(!issueConv)]
  print("Done with fitting terminal models")
  print(paste0("Proportion of runs removed due to issues: ", mean(issueConv)))


  attributes(fitsim)$simN = simN
  attributes(fitsim)$simF = simF
  attributes(fitsim)$fit = fit

  return(fitsim)
}

nsim = 4000 #Note, may run with lower nsim and combine several runs later
ncores = 10
fitsim = simFullData(fitOriginal,ncores = ncores,
                     seed = 12345,nsim = nsim,
                     rec.years = fitOriginal$data$years)
rets = lapply(1:length(fitsim),function(i){
  out <- tryCatch(
    {
      print(i)
      retro(fitsim[[i]],year = 5,ncores = 5);
    },
    error=function(cond) {
      return(NA)
    }
  )
  out
})

issueRet = rep(FALSE,length(rets))
for(i in 1:length(rets)){
  if(length(rets[[i]])==1){
    issueRet[i] = TRUE
  }else{
    for(ii in 1:length(rets[[i]])){
      if(rets[[i]][[ii]]$opt$convergence != 0){
        issueRet[i] = TRUE
      }
    }
  }
}

print("Done with fitting terminal models")
print(paste0("Proportion of runs removed due to issues: ", round(mean(issueRet),3)))

rets = rets[which(!issueRet)]



fit = attributes(rets[[1]])$fit
lagF = max(fit$data$years) -max(as.numeric(rownames(getFleet(fit,1))))

if(is.null(nsim)) nsim = length(rets)
mohR = mohn(rets[[1]])[1]
mohSSB = mohn(rets[[1]])[2]
mohF = mohn(rets[[1]],lag = lagF)[3]
moh = c(mohR,mohSSB,mohF)
for(i in 2:length(rets)){
  mohR = mohn(rets[[i]])[1]
  mohSSB = mohn(rets[[i]])[2]
  mohF = mohn(rets[[i]],lag = lagF)[3]
  mohTmp = c(mohR,mohSSB,mohF)
  moh = rbind(moh,mohTmp)
}
moh = as.data.frame(moh)
names(moh) = c("R", "SSB","FBar")

ssb = rep(0,length(rets))
fbar = rep(0,length(rets))
sdN = rep(0,length(rets))
sdF = rep(0,length(rets))
for(i in 1:length(rets)){
  ssb[i] = rev(ssbtable(rets[[i]][[5]])[,1])[1]
  fbar[i] = rev(fbartable(rets[[i]][[5]])[,1])[1]
  rlSd = as.list(rets[[i]][[5]]$sdrep,what = "Std",report = TRUE)
  sdN[i] = rev(rlSd$logssb)[1]
  sdF[i] = rev(rlSd$logfbar)[2]
}

moh$ssb = ssb
moh$fbar = fbar
moh$sdN = sdN
moh$sdF = sdF



if(FALSE){#Save samples with assosiated variables
  mohTOld = read.table("results/NEACodFullSimDetailed.txt", header = TRUE)
  moh = rbind(mohTOld,moh)
  write.table(moh,file = "results/NEACodFullSimDetailed.txt", row.names = FALSE)
}


moh = read.table("results/NEACodFullSimDetailed.txt", header = TRUE)[1:nsim,]

plot(moh$ssb,moh$SSB)
plot(moh$fbar,moh$SSB)
plot(moh$ssb,moh$FBar)
plot(moh$fbar,moh$FBar)
plot(moh$sdN,moh$SSB,xlim = c(0.05,0.20))
plot(moh$sdF,moh$FBar)

plot(moh$ssb,moh$sdN)

nL = which(moh$sdN<= sort(moh$sdN)[1000])
nU = which(moh$sdN>= sort(moh$sdN)[dim(moh)[1]-1001])
quantile(moh$SSB[nL], c(0.025,0.975))
quantile(moh$SSB[nU], c(0.025,0.975))

fL = which(moh$sdF<= sort(moh$sdF)[1000])
fU = which(moh$sdF>= sort(moh$sdF)[dim(moh)[1]-1001])
quantile(moh$FBar[fL], c(0.025,0.975))
quantile(moh$FBar[fU], c(0.025,0.975))


ssbL = which(moh$ssb< quantile(moh$ssb,0.25, na.rm = TRUE))
ssbU = which(moh$ssb> quantile(moh$ssb,0.75, na.rm = TRUE))
quantile(moh$SSB[ssbL], c(0.025,0.975))
quantile(moh$SSB[ssbU], c(0.025,0.975))

fbarL = which(fbar< quantile(fbar,0.25, na.rm = TRUE))
fbarU = which(fbar> quantile(fbar,0.75, na.rm = TRUE))
quantile(mohT$SSB[fbarL], c(0.025,0.975))
quantile(mohT$SSB[fbarU], c(0.025,0.975))



#Original uncertainty
rlOrig = as.list(retOriginal[[5]]$sdrep,what = "Std",report = TRUE)
rev(rlOrig$logfbar)[2]
rev(rlOrig$logssb)[1]



###########################################

#Illustrate that we are able to estimate the model with this simulation procedure also.
nn = length(fitOriginal$data$years)

x11(width =15,height = 15)
par(mfrow = c(4,4))
for(ii in 1:16){
  tmp = (dim(attributes(fitsim)$fit$pl$logN)[2]-nn):dim(attributes(fitsim)$fit$pl$logN)[2]
  ssb = calcSSB(attributes(fitsim)$simN[[ii]], attributes(fitsim)$fit$data$stockMeanWeight[tmp,], attributes(fitsim)$fit$data$propMat[tmp,],
                attributes(fitsim)$fit$data$natMor[tmp,], attributes(fitsim)$fit$data$propF[tmp,],
                attributes(fitsim)$fit$data$propM[tmp,], attributes(fitsim)$simF[[ii]])

  ssbplot(fitsim[[ii]],main = paste0("SSB, simulation ",ii))
  points(attributes(fitsim)$fit$data$years,ssb,type = 'l', col = 2)

}


x11(width =15,height = 15)
par(mfrow = c(4,4))
for(ii in 1:16){
  fbar = calcFbar(attributes(fitsim)$simF[[ii]], ageRange = attributes(fitsim)$fit$conf$fbarRange)

  fbarplot(fitsim[[ii]],main = paste0("Fbar, simulation ",ii))
  points(attributes(fitsim)$fit$data$years,fbar,type = 'l', col = 2)

}



