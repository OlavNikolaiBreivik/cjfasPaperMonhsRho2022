library(parallel)


addrow<-function(mat,v){
  rn <- rownames(mat)
  mat <-rbind(mat,v)
  rownames(mat)<-c(rn,max(as.integer(rn))+1)
  return(mat)
}

calcFbar = function(F, ageRange){
  FF = F[which(as.numeric(rownames(F)) %in% ageRange[1]:ageRange[2]),]
  return(colMeans(FF))
}

calcSSB = function(N,sw = sw, mo = mo, nm = NULL,pf = NULL,pm = NULL,F = NULL){
  if(!is.null(pf)){
    zz = exp(-F*t(pf) -t(nm)*t(pm))
    ssb = N*t(sw)*t(mo)*zz#TODO: pf and pm
  }else{
    ssb = N*t(sw)*t(mo)#TODO: pf and pm
  }
  return(colSums(ssb))
}

#processNoiseF = TRUE;ncores = 1;  scaleQ = rep(1,length(ret)+1);scaleC= rep(1,length(ret)+1);scaleM= rep(1,length(ret)+1);scaleSdF= rep(1,length(ret)+1);scaleSdObsC= rep(1,length(ret)+1);scaleSdObsI= rep(1,length(ret)+1);scaleCorObs= FALSE;corF = FALSE;whichQScale = seq(1,999)
mohnsDist = function(ret,nsim = 50,seed = NULL, processNoiseF = TRUE,ncores = parallel::detectCores(),
                     scaleQ = rep(1,length(ret)+1),scaleC= rep(1,length(ret)+1),scaleM= rep(1,length(ret)+1),scaleSdF= rep(1,length(ret)+1),scaleSdObsC= rep(1,length(ret)+1),
                     scaleSdObsI= rep(1,length(ret)+1),scaleCorObs= FALSE,corF = FALSE,whichQScale = seq(1,999),...){
  nya = length(ret)
  if(!is.null(seed)){
    set.seed(seed)
  }

  simDataList<-list()
  simN = list()
  simF = list()

  fit = ret[[length(ret)]]
  fitTerminal = attributes(ret)$fit
  nm <- fitTerminal$data$natMor
  cw <- fitTerminal$data$catchMeanWeight
  dw <- fitTerminal$data$disMeanWeight
  lw <- fitTerminal$data$landMeanWeight
  lf <- fitTerminal$data$landFrac
  mo <- fitTerminal$data$propMat
  sw <- fitTerminal$data$stockMeanWeight
  pf <- fitTerminal$data$propF
  pm <- fitTerminal$data$propM
  nn = dim(nm)[2]

  fit$data$natMor = fitTerminal$data$natMor
  fc <- forecast(fit, fscale=rep(1,nya+1), nosim=nsim, processNoiseF=processNoiseF, savesim = TRUE,
                 scaleM = scaleM,scaleSdF = scaleSdF,corF = corF,...)



  for(j in 1:nsim){
    cnnew <- getFleet(fit,1)
    surveysnew <- list()
    for(ff in 2:fit$data$noFleets){
      surveysnew[[ff-1]] <- getFleet(fit,ff)
      attr(surveysnew[[ff-1]],"time") <- rep(fit$data$sampleTimes[ff],2)
    }

    lag = rep(0, fit$data$noFleets)
    for(ff in 1:fit$data$noFleets){
      lag[ff] = max(fit$data$years) -max(as.numeric(rownames(getFleet(fit,ff))))
    }

    simN[[j]] = matrix(0,nrow = nn, ncol = nya+1)
    simF[[j]] = matrix(0,nrow = nn, ncol = nya+1)
    rownames(simF[[j]]) = fit$conf$minAge:fit$conf$maxAge
    rownames(simN[[j]]) = fit$conf$minAge:fit$conf$maxAge
    simN[[j]][,1] = exp(fc[[1]]$sim[j,1:nn])
    simF[[j]][,1] = exp(fc[[1]]$sim[j,attr(fc,"fit")$conf$keyLogFsta[1,]+(nn+1)])
    for(y in 2:(nya+1)){
      nmThis = nm[length(fit$data$years)-1 + y,]*scaleM[y]
      nmPrev = nm[length(fit$data$years)-1 + y-1,]*scaleM[y-1]

      swThis = sw[length(fit$data$years)-1 + y,]
      swPrev = sw[length(fit$data$years)-1 + y-1,]
      moThis = mo[length(fit$data$years)-1 + y,]
      moPrev = mo[length(fit$data$years)-1 + y-1,]

      theFleets<-list()
      N<-exp(fc[[y]]$sim[j,1:nn])
      F<-exp(fc[[y]]$sim[j,attr(fc,"fit")$conf$keyLogFsta[1,]+(nn+1)])
      Nprev<-exp(fc[[y-1]]$sim[j,1:nn])
      Fprev<-exp(fc[[y-1]]$sim[j,attr(fc,"fit")$conf$keyLogFsta[1,]+(nn+1)])

      simN[[j]][,y] = N
      simF[[j]][,y] = F
      if(lag[1]==0){
        Cpred <- F/(F+nmThis)*(1-exp(-F-nmThis))*N
      }else if(lag[1]==1){
        Cpred <- Fprev/(Fprev+nmPrev)*(1-exp(-Fprev-nmPrev))*Nprev
      }else{
        stop("Lag larger than 1 for catch")
      }

      covCatch <- attr(fc,"fit")$rep$obsCov[[1]]
      if(lag[1]==0){
        Csim <- exp(log(Cpred)+rmvnorm(1,rep(0,nrow(covCatch)),covCatch*scaleSdObsC[y]^2))*scaleC[y]
      }else if(lag[1]==1){
        Csim <- exp(log(Cpred)+rmvnorm(1,rep(0,nrow(covCatch)),covCatch*scaleSdObsC[y-1]^2))*scaleC[y-1]
      }

      theFleets[[1]] <- Csim

      f <- Vectorize(function(ii)ifelse(ii==-1,NA,attr(fc,"fit")$opt$par[ii+1]))
      Q <- t(exp(apply(attr(fc,"fit")$conf$keyLogFpar,1,f)))

      Qtmp = Q #Applied when defining plus group.
      for(ff in 2:fit$data$noFleets){
        if(fit$conf$maxAgePlusGroup[ff]==1){
          Q[ff,which(is.na(Q[ff,]))] = Q[ff,max(which(Q[ff,]>0))]
        }
        if(ff %in% whichQScale){
          if(lag[ff]==0){
            Q[ff,] <- Q[ff,]*scaleQ[y]
          }else if(lag[ff]==1){
            Q[ff,] <- Q[ff,]*scaleQ[y-1]
          }
        }

      }

      tau <- attr(fc,"fit")$data$sampleTimes

      for(ff in 2:fit$data$noFleets){
        Qf <- Q[ff,]
        if(fit$data$fleetTypes[[ff]]==2){
          Sfpred <- Qf*N*exp(-tau[ff]*(F+nmThis))
          SfpredPrev <- Qf*Nprev*exp(-tau[ff]*(Fprev+nmPrev))
          if(fit$conf$maxAgePlusGroup[ff]==1){
            Sfpred[max(which(Qtmp[ff,]>0))] = sum(Sfpred[max(which(Qtmp[ff,]>0)): nn])
            SfpredPrev[max(which(Qtmp[ff,]>0))] = sum(SfpredPrev[max(which(Qtmp[ff,]>0)): nn])
          }

          if(sum(is.na(Qtmp[ff,]))>0){
            Sfpred <- Sfpred[-which(is.na(Qtmp[ff,]))]
            SfpredPrev <- SfpredPrev[-which(is.na(Qtmp[ff,]))]
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
            SfpredPrev <- Qf[1]*sum(N*exp(-tau[ff]*(F+nmThis))*swPrev*moPrev)
          }else if(fit$conf$keyBiomassTreat[[ff]]==5){
            #TSB-index #NB!!! Can't see propM og propF being used (but applied for SSB predictions)
            Sfpred <- Qf[1]*sum(N*exp(-tau[ff]*(F+nmThis))*swThis)
            SfpredPrev <- Qf[1]*sum(N*exp(-tau[ff]*(F+nmPrev))*swPrev)
          }
          else if(fit$conf$keyBiomassTreat[[ff]]==6){
            #TSN-index
            Sfpred <- Qf[1]*sum(N*exp(-tau[ff]*(F+nmThis)))
            SfpredPrev <- Qf[1]*sum(N*exp(-tau[ff]*(F+nmThis)))
          }else{
            stop(paste0("Mohns dist not yet implemented for bimoassTreat ",fit$conf$keyBiomassTreat[[ff]]))
          }
          Sfcov = as.matrix(exp(fit$pl$logSdLogObs[fit$conf$keyVarObs[ff,1]+1])^2)
        }
        if(lag[ff]==0){
          Sfsim <- exp(log(Sfpred)+rmvnorm(1,rep(0,nrow(Sfcov)),Sfcov*scaleSdObsI[y]^2))
        }else if(lag[ff]==1){
          Sfsim <- exp(log(SfpredPrev)+rmvnorm(1,rep(0,nrow(Sfcov)),Sfcov*scaleSdObsI[y-1]^2))
        }else{
          Sfsim = rep(NA,nrow(Sfcov))
          warning("Lag larger than 1, set to NA")
        }
        theFleets[[ff]] <- Sfsim
      }

      cnnew <- addrow(cnnew,theFleets[[1]])
      for(ff in 2:fit$data$noFleets){
        surveysnew[[ff-1]] <- addrow(surveysnew[[ff-1]],theFleets[[ff]])
      }
    }
    for(ff in 2:fit$data$noFleets){

      ss = surveysnew[[ff-1]]
      ssOrig = getFleet(fitTerminal,ff)

      if(min(dim(ss)==dim(ssOrig))==1){
        ss[which(is.na(ssOrig))] = NA
      }else{
        stop("Did not end simulations at terminal year")
      }
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

  rets= lapply(1:length(fitsim), function(i) {
    print(i)
    out = tryCatch({
        retro(fitsim[[i]],year = nya,ignore.parm.uncertainty=TRUE, newtonsteps = 0,ncores = nya)
      },
      error=function(cond){
        return(NA)
      })
    out
    }
  )

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

  simN = simN[which(!issueRet)]
  simF = simF[which(!issueRet)]
  rets = rets[which(!issueRet)]

  attributes(rets)$simN = simN
  attributes(rets)$simF = simF
  attributes(rets)$fit = fitTerminal
  attributes(rets)$fitsim = fitsim[which(!issueRet)]

  class(rets) = "mohnsDist"

  return(rets)
}

plotMohnDist = function(rets,ret=NULL, what = "mohn", examples = FALSE,mohT = NULL,p = 0.05,main ="",...){

  fit = attributes(rets[[1]])$fit
  lagF = max(fit$data$years) -max(as.numeric(rownames(getFleet(fit,1))))

  if(lagF>1){
    stop("Catch not given in terminal year or terminal year-1")
  }

  nya = length(rets[[1]])
  nsim = length(rets)
  if(!examples){
    mohR = mohn(rets[[1]])[1]
    mohSSB = mohn(rets[[1]])[2]
    mohF = mohn(rets[[1]],lag = lagF)[3]
    moh = c(mohSSB,mohF,mohR)
    for(i in 2:nsim){
      mohR = mohn(rets[[i]])[1]
      mohSSB = mohn(rets[[i]])[2]
      mohF = mohn(rets[[i]],lag = lagF)[3]
      mohTmp = c(mohSSB,mohF,mohR)
      moh = rbind(moh,mohTmp)
    }

    qR = quantile(moh[,1], c(p/2,1-p/2))
    qSSB = quantile(moh[,2], c(p/2,1-p/2))
    qF = quantile(moh[,3], c(p/2,1-p/2))

    ylim = c(min(moh), max(moh))
    library(plotrix)
    plotCI(1,y = median(moh[,1]), li = qR[1], ui = qR[2], ylim = ylim, xlim = c(0.5,3.3),
           main = main,
           xlab = "",
           ylab = "Mohn's rho",
           xaxt='n', lwd = 2, cex.main = 2)
    for(i in 2:dim(moh)[1]){
      if(moh[i,1] > qR[2] | moh[i,1] < qR[1]) points(1,moh[i,1], cex = 0.5)
    }

    plotCI(2,y = median(moh[,2]), li = qSSB[1], ui = qSSB[2], ylim = ylim, lwd = 2, add = TRUE)
    for(i in 2:dim(moh)[1]){
      if(moh[i,2] > qSSB[2] | moh[i,2] < qSSB[1]) points(2,moh[i,2], cex = 0.5)
    }

    plotCI(3,y = median(moh[,3]), li = qF[1], ui = qF[2], ylim = ylim, lwd = 2, add = TRUE)
    for(i in 2:dim(moh)[1]){
      if(moh[i,3] > qF[2] | moh[i,3] < qF[1]) points(3,moh[i,3], cex = 0.5)
    }

    abline(h = 0)

    axis(1, at=1:3, labels=c("SSB","F-bar","R"), las=1)
    mohRObs = mohn(ret)[1]
    mohSSBObs = mohn(ret)[2]
    mohFObs = mohn(ret,lag = lagF)[3]
    mohObs = c(mohSSBObs,mohFObs,mohRObs)
    points(mohObs, col = 'red', cex = 1, pch = 16)


    pp = rep(0,3)
    for(i in 1:3){
      pp[i] = round(min(mean(moh[,i]>mohObs[i]),mean(moh[,i]<mohObs[i]))*2,2)
    }
    if(pp[1]<p/2) col = 'red' else col = "black"
    text(0.7, 0.02, paste0("P:  ",pp[1]),
         cex=1.3, pos=3,col=col)
    if(pp[2]<p/2) col = 'red' else col = "black"
    text(1.7, 0.02, paste0("P:  ",pp[2]),
         cex=1.3, pos=3,col=col)
    if(pp[3]<p/2) col = 'red' else col = "black"
    text(2.7, 0.02, paste0("P:  ",pp[3]),
         cex=1.3, pos=3,col=col)

    points(c(0.6,1.4),c(0.2,0.2), type = 'l', lty = 2,col = 'green',lwd = 2)
    points(c(0.6,1.4),c(-0.15,-0.15), type = 'l', lty = 2,col = 'green',lwd = 2)
  }

  if(examples){
    if(is.null(ret))stop("need ret as input")
    par(mfrow = c(1,3))

    fitFirst = ret[[length(ret)]]
    fbar = calcFbar(attributes(rets)$simF[[1]], ageRange = attributes(rets)$fit$conf$fbarRange)
    tmp = (dim(attributes(rets)$fit$pl$logN)[2]-nya):dim(attributes(rets)$fit$pl$logN)[2]
    ssb = calcSSB(attributes(rets)$simN[[1]], attributes(rets)$fit$data$stockMeanWeight[tmp,], attributes(rets)$fit$data$propMat[tmp,],
                  attributes(rets)$fit$data$natMor[tmp,], attributes(rets)$fit$data$propF[tmp,],
                  attributes(rets)$fit$data$propM[tmp,], attributes(rets)$simF[[1]])

    ff = attributes(rets[[1]])$fit
    ylimF = c(0, max(fbartable(ff)[(dim(fbartable(ff))[1]-10):dim(fbartable(ff))[1],3]))
    ylimSSB = c(0, max(ssbtable(ff)[(dim(ssbtable(ff))[1]-10):dim(ssbtable(ff))[1],3]))
    ylimR = c(0, max(rectable(ff)[(dim(rectable(ff))[1]-10):dim(rectable(ff))[1],3]))

#    ssbplot(rets[[1]],main = paste0("SSB, simulation ",1),ylim = ylimSSB,...)
    ssbplot(rets[[1]],main = paste0("SSB, simulation ",1),...)
    points(max(fitFirst$data$years) + 0:(length(tmp)-1),ssb,type = 'l', col = 2)
    oldpar = par(ask = TRUE); on.exit(par(oldpar))

    fit = attributes(ret)$fit
    drop=max(fit$data$aux[,"year"])-max(fit$data$aux[fit$data$aux[,"fleet"]==1,"year"])
#    fbarplot(rets[[1]],main = paste0("Fbar, simulation ",1),ylim = ylimF,...)
    fbarplot(rets[[1]],main = paste0("Fbar, simulation ",1),...)
    if(drop ==0){
      points(max(fitFirst$data$years) + 0:(length(tmp)-1),fbar,type = 'l', col = 2)
    }else if(drop==1){
      points(max(fitFirst$data$years) + 0:(length(tmp)-2),fbar[-length(fbar)],type = 'l', col = 2)
    }else{
      stop("Catch not given in terminal year or terminal year -1")
    }

#    recplot(rets[[1]],main = paste0("R, simulation ",1),ylim = ylimR,...)
    recplot(rets[[1]],main = paste0("R, simulation ",1),...)
    points(max(fitFirst$data$years) + 0:(length(tmp)-1),attributes(rets)$simN[[1]][1,],type = 'l', col = 2)

    for(ii in 2:length(rets)){
      fbar = calcFbar(attributes(rets)$simF[[ii]], ageRange = attributes(rets)$fit$conf$fbarRange)
      tmp = (dim(attributes(rets)$fit$pl$logN)[2]-nya):dim(attributes(rets)$fit$pl$logN)[2]
      ssb = calcSSB(attributes(rets)$simN[[ii]], attributes(rets)$fit$data$stockMeanWeight[tmp,], attributes(rets)$fit$data$propMat[tmp,],
                    attributes(rets)$fit$data$natMor[tmp,], attributes(rets)$fit$data$propF[tmp,],
                    attributes(rets)$fit$data$propM[tmp,], attributes(rets)$simF[[ii]])

      ff = attributes(rets[[ii]])$fit
      ylimF = c(0, max(fbartable(ff)[(dim(fbartable(ff))[1]-10):dim(fbartable(ff))[1],3]))
      ylimSSB = c(0, max(ssbtable(ff)[(dim(ssbtable(ff))[1]-10):dim(ssbtable(ff))[1],3]))
      ylimR = c(0, max(rectable(ff)[(dim(rectable(ff))[1]-10):dim(rectable(ff))[1],3]))

#      ssbplot(rets[[ii]],main = paste0("SSB, simulation ",ii),ylim = ylimSSB,...)
      ssbplot(rets[[ii]],main = paste0("SSB, simulation ",ii),...)
      points(max(fitFirst$data$years) + 0:(length(tmp)-1),ssb,type = 'l', col = 2)
#      fbarplot(rets[[ii]],main = paste0("Fbar, simulation ",ii),ylim = ylimF,...)
      fbarplot(rets[[ii]],main = paste0("Fbar, simulation ",ii),...)
      if(drop ==0){
        points(max(fitFirst$data$years) + 0:(length(tmp)-1),fbar,type = 'l', col = 2)
      }else if(drop==1){
        points(max(fitFirst$data$years) + 0:(length(tmp)-2),fbar[-length(fbar)],type = 'l', col = 2)
      }
#      recplot(rets[[ii]],main = paste0("R, simulation ",ii),ylim =  ylimR, ...)
      recplot(rets[[ii]],main = paste0("R, simulation ",ii), ...)
      points(max(fitFirst$data$years) + 0:(length(tmp)-1),attributes(rets)$simN[[ii]][1,],type = 'l', col = 2)
    }
  }
}


saveMohnsSamples = function(rets,file,nsim = NULL){

  fit = attributes(rets[[1]])$fit
  lagF = max(fit$data$years) -max(as.numeric(rownames(getFleet(fit,1))))

  if(is.null(nsim)) nsim = length(rets)
  mohR = mohn(rets[[1]])[1]
  mohSSB = mohn(rets[[1]])[2]
  mohF = mohn(rets[[1]],lag = lagF)[3]
  moh = c(mohR,mohSSB,mohF)
  for(i in 2:nsim){
    mohR = mohn(rets[[i]])[1]
    mohSSB = mohn(rets[[i]])[2]
    mohF = mohn(rets[[i]],lag = lagF)[3]
    mohTmp = c(mohR,mohSSB,mohF)
    moh = rbind(moh,mohTmp)
  }
  write.table(moh,file = file, row.names = FALSE)
}






#t-1 in getAve for natmor and include model violation options
forecast <- function(fit, fscale=NULL, catchval=NULL, catchval.exact=NULL, fval=NULL, nextssb=NULL, landval=NULL, cwF=NULL, nosim=1000, year.base=max(fit$data$years), ave.years=max(fit$data$years)+(-4:0), rec.years=max(fit$data$years)+(-9:0), label=NULL, overwriteSelYears=NULL, deterministic=FALSE, processNoiseF=TRUE, customWeights=NULL, customSel=NULL, lagR=FALSE, splitLD=FALSE, addTSB=FALSE, useSWmodel=(fit$conf$stockWeightModel==1), useCWmodel=(fit$conf$catchWeightModel==1), useMOmodel=(fit$conf$matureModel==1), useNMmodel=(fit$conf$mortalityModel==1), savesim=FALSE,
                     scaleM= rep(1,999),scaleSdF= rep(1,999),corF = FALSE,recPool= NULL,scalInitialNsd = 1){

  resample <- function(x, ...){
    if(deterministic){
      ret <- mean(x)
    }else{
      ret <- x[sample.int(length(x), ...)]
    }
    return(ret)
  }
  ns<-max(length(fscale), length(catchval), length(catchval.exact), length(fval), length(nextssb), length(cwF))
  if(missing(fscale)&missing(fval)&missing(catchval)&missing(catchval.exact)&missing(nextssb)&missing(cwF))stop("No scenario is specified")
  if(missing(fscale)) fscale <- rep(NA,ns)
  if(missing(fval)) fval <- rep(NA,ns)
  if(missing(catchval)) catchval <- rep(NA,ns)
  if(missing(catchval.exact)) catchval.exact <- rep(NA,ns)
  if(missing(nextssb)) nextssb <-rep(NA,ns)
  if(missing(landval)) landval <-rep(NA,ns)
  if(missing(cwF)) cwF <-rep(NA,ns)

  if(!all(rowSums(!is.na(cbind(fscale, catchval, catchval.exact, fval, nextssb, landval, cwF)))==1)){
    stop("For each forecast year exactly one of fscale, catchval or fval must be specified (all others must be set to NA)")
  }

  if(!is.null(overwriteSelYears)){
    fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)
    Ftab <- faytable(fit)
    fixedsel <- colMeans(Ftab[as.integer(rownames(Ftab))%in%overwriteSelYears,,drop=FALSE])
    fixedsel <- fixedsel/mean(fixedsel[fromto[1]:fromto[2]])
  }

  if(!is.null(customSel)){
    fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)
    customSel <- customSel/mean(customSel[fromto[1]:fromto[2]])
  }

  odat<-fit$obj$env$data
  opar<-fit$obj$env$parList(par = fit$obj$env$last.par.best)
  omap<-fit$obj$env$map
  omap[names(which(lapply(opar,length)==0))]<-NULL
  if(useSWmodel & (fit$conf$stockWeightModel==0)){
    stop("stockWeightModel cannot be used for forecast when it was not part of the fitted model")
  }
  if(useSWmodel){
    rnSW<-1:(nrow(opar$logSW)+ns)+as.integer(rownames(odat$stockMeanWeight)[1])-1
    opar$logSW<-matrix(0,nrow=nrow(opar$logSW)+ns, ncol=ncol(opar$logSW))
    oran<-unique(names(fit$obj$env$par[fit$obj$env$random]))
    objsw <- MakeADFun(odat, opar, random = oran, DLL = "stockassessment", map=omap)
    sdrep<- sdreport(objsw, par.fixed=fit$opt$par, ignore.parm.uncertainty=TRUE)
    idx<-names(sdrep$value)=="logSW"
    simLogSw <-rmvnorm(nosim, sdrep$value[idx], sdrep$cov[idx,idx])
  }
  if(useCWmodel & (fit$conf$catchWeightModel==0)){
    stop("catchWeightModel cannot be used for forecast when it was not part of the fitted model")
  }
  if(useCWmodel){
    rnCW<-1:(nrow(opar$logCW)+ns)+as.integer(rownames(odat$catchMeanWeight)[1])-1
    opar$logCW<-matrix(0,nrow=nrow(opar$logCW)+ns, ncol=ncol(opar$logCW))
    oran<-unique(names(fit$obj$env$par[fit$obj$env$random]))
    objcw <- MakeADFun(odat, opar, random = oran, DLL = "stockassessment", map=omap)
    sdrep<- sdreport(objcw, par.fixed=fit$opt$par, ignore.parm.uncertainty=TRUE)
    idx<-names(sdrep$value)=="logCW"
    simLogCW <-rmvnorm(nosim, sdrep$value[idx], sdrep$cov[idx,idx])
  }
  if(useNMmodel & (fit$conf$mortalityModel==0)){
    stop("mortalityModel cannot be used for forecast when it was not part of the fitted model")
  }
  if(useNMmodel){
    rnNM<-1:(nrow(opar$logNM)+ns)+as.integer(rownames(odat$natMor)[1])-1
    opar$logNM<-matrix(0,nrow=nrow(opar$logNM)+ns, ncol=ncol(opar$logNM))
    oran<-unique(names(fit$obj$env$par[fit$obj$env$random]))
    objnm <- MakeADFun(odat, opar, random = oran, DLL = "stockassessment", map=omap)
    sdrep<- sdreport(objnm, par.fixed=fit$opt$par, ignore.parm.uncertainty=TRUE)
    idx<-names(sdrep$value)=="logNM"
    simLogNM <-rmvnorm(nosim, sdrep$value[idx], sdrep$cov[idx,idx])
  }
  if(useMOmodel & (fit$conf$matureModel==0)){
    stop("matureModel cannot be used for forecast when it was not part of the fitted model")
  }
  if(useMOmodel){
    rnMO<-1:(nrow(opar$logitMO)+ns)+as.integer(rownames(odat$propMat)[1])-1
    opar$logitMO<-matrix(0,nrow=nrow(opar$logitMO)+ns, ncol=ncol(opar$logitMO))
    oran<-unique(names(fit$obj$env$par[fit$obj$env$random]))
    objmo <- MakeADFun(odat, opar, random = oran, DLL = "stockassessment", map=omap)
    sdrep<- sdreport(objmo, par.fixed=fit$opt$par, ignore.parm.uncertainty=TRUE)
    idx<-names(sdrep$value)=="logitMO"
    simLogitMO <-rmvnorm(nosim, sdrep$value[idx], sdrep$cov[idx,idx])
  }

  getF <- function(x, allowSelOverwrite=FALSE){
    idx <- fit$conf$keyLogFsta[1,]+1
    nsize <- length(idx)
    ret <- exp(x[nsize+idx])
    ret[idx==0] <- 0
    if(allowSelOverwrite){
      if(!is.null(overwriteSelYears)){
        fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)
        thisfbar<-mean(ret[fromto[1]:fromto[2]])
        ret<-fixedsel*thisfbar
      }
      if(!is.null(customSel)){
        fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)
        thisfbar<-mean(ret[fromto[1]:fromto[2]])
        ret<-customSel*thisfbar
      }
    }
    ret
  }

  fbar<-function(x){
    fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)
    mean(getF(x)[fromto[1]:fromto[2]])
  }

  fbarFrac<-function(x, lf){
    fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)
    mean((lf*getF(x))[fromto[1]:fromto[2]])
  }

  getCWF<-function(x,w){
    sum(getF(x)*w)
  }

  getN <- function(x){
    idx <- fit$conf$keyLogFsta[1,]+1
    nsize <- length(idx)
    ret <- exp(x[1:nsize])
    ret
  }

  getState <- function(N,F){
    k <- fit$conf$keyLogFsta[1,]
    F <- F[k>=0]
    k <- !duplicated(k[k>=0])
    x <- log(c(N,F[k]))
    x
  }

  getProcessVar <- function(fit,i){
    cof <- coef(fit)
    sdN <- exp(cof[names(cof)=="logSdLogN"][fit$conf$keyVarLogN+1])
    sdN[1]<-0
    nN <- length(sdN)
    sdF <- exp(cof[names(cof)=="logSdLogFsta"][fit$conf$keyVarF+1])*scaleSdF[i+1]
    if(processNoiseF==FALSE){
      sdF <- sdF*0
    }
    k<-fit$conf$keyLogFsta[1,]
    sdF <- sdF[k>=0]
    k <- unique(k[k >= 0] + 1)
    sdF <-sdF[k]
    nF <- length(sdF)
    if(fit$conf$corFlag==0){
      corr <- diag(nF)
    }
    if(fit$conf$corFlag==1){
      y <- cof[names(cof)=="itrans_rho"]
      rho <- 2/(1+exp(-2*y))-1
      corr <- matrix(rho, nrow=nF, ncol=nF)
      diag(corr) <- 1
    }
    if(fit$conf$corFlag==2){
      y <- cof[names(cof)=="itrans_rho"]
      rho <- 2/(1+exp(-2*y))-1
      corr <- diag(sdF)
      corr <- rho^abs(row(corr)-col(corr))
    }
    if(corF){
      corr <- diag(nF)
      if(fit$conf$corFlag==0)warning("Try to switch of corr F which is already off")
    }
    cov <- matrix(0,nrow=nN+nF,ncol=nN+nF)
    cov[1:nN,1:nN] <- diag(sdN^2)

    if(fit$conf$corFlag <3){
      cov[nN+1:nF,nN+1:nF] <- (sdF%*%t(sdF))*corr
    }else{
      if(fit$conf$corFlag ==3){
        sdU <- exp(cof[names(cof)=="sepFlogSd"][1])
        sdV <- exp(cof[names(cof)=="sepFlogSd"][2])

        diag(cov[nN+1:(nF-1),nN+1:(nF-1)]) <- sdU^2
        cov[nN+nF,nN+nF] <- sdV^2
      }
    }
    cov
  }

  step <- function(x, nm, recpool, scale, inyear=FALSE){
    F <- getF(x, allowSelOverwrite=!inyear)
    N <- getN(x)
    if(!inyear){
      Z <- F+nm
      n <- length(N)
      N <- c(resample(recpool,1),N[-n]*exp(-Z[-n])+c(rep(0,n-2),N[n]*exp(-Z[n])))
    }
    F <- F*scale
    xx <- getState(N,F)
    return(xx)
  }

  scaleF <- function(x, scale){
    F <- getF(x)*scale
    N <- getN(x)
    xx <- getState(N,F)
    return(xx)
  }

  sel<-function(x){
    getF(x)/fbar(x)
  }

  catch <- function(x, nm, cw){
    F <- getF(x)
    Z <- F+nm
    N <- getN(x)
    C <- F/Z*(1-exp(-Z))*N
    return(sum(cw*C))
  }

  catchFrac <- function(x, nm, w, frac){
    F <- getF(x)
    Z <- F+nm
    N <- getN(x)
    C <- F/Z*(1-exp(-Z))*N
    return(sum(frac*w*C))
  }

  catchatage <- function(x, nm){
    F <- getF(x)
    Z <- F+nm
    N <- getN(x)
    C <- F/Z*(1-exp(-Z))*N
    return(C)
  }

  ssb <- function(x, nm, sw, mo, pm, pf){
    F <- getF(x)
    N <- getN(x)*exp(-pm*nm-pf*F)
    return(sum(N*mo*sw))
  }

  tsb <- function(x, sw){
    F <- getF(x)
    N <- getN(x)
    return(sum(N*sw))
  }

  rectab<-rectable(fit)

  if(missing(recPool)){
    recpool<-rectab[rownames(rectab)%in%rec.years,1]
  }else{
    recpool = recPool#Used when comparing Mohn's rho from different scenarios
  }

  # Get final state
  if(year.base==max(fit$data$years)){
    est <- fit$sdrep$estY
    cov <- fit$sdrep$covY
  } else if(year.base==(max(fit$data$years)-1)){
    est <- fit$sdrep$estYm1
    cov <- fit$sdrep$covYm1
  }else if(year.base == min(fit$data$years)){
    est <- fit$sdrep$estYfirst
    cov <- fit$sdrep$covYfirst
  }else{
    stop("State not saved, so cannot proceed from this year")
  }

  if(scalInitialNsd !=1){
    M = diag(dim(cov)[1])
    idx <- fit$conf$keyLogFsta[1,]+1
    nsize <- length(idx)
    M[1:nsize,1:nsize] = M[1:nsize,1:nsize] *scalInitialNsd
    cov = M%*%cov%*%M
  }

  if(deterministic)cov<-cov*0
  sim<-rmvnorm(nosim, mu=est, Sigma=cov)

  if(is.null(overwriteSelYears) & is.null(customSel)){
    if(!isTRUE(all.equal(est,getState(getN(est),getF(est)))))stop("Sorry somthing is wrong here (check code for getN, getF, and getState)")
  }

  doAve <- function(x,y)colMeans(x[rownames(x)%in%ave.years,,drop=FALSE], na.rm=TRUE)
  ave.sw <- doAve(fit$data$stockMeanWeight)
  ave.cw <- doAve(fit$data$catchMeanWeight)
  ave.mo <- doAve(fit$data$propMat)
  ave.nm <- doAve(fit$data$natMor)
  ave.lf <- doAve(fit$data$landFrac)
  ave.lw <- doAve(fit$data$landMeanWeight)
  ave.pm <- doAve(fit$data$propM)
  ave.pf <- doAve(fit$data$propF)
  getThisOrAve <- function(x,y, ave){
    if((y %in% rownames(x))&(!all(is.na(x[which(rownames(x)==y),])))){
      ret <- x[which(rownames(x)==y),]
    }else{
      ret <- ave
    }
    ret
  }
#  procVar<-getProcessVar(fit)
  simlist<-list()
  rs <- NULL
  for(i in 0:(length(fscale)-1)){
    procVar<-getProcessVar(fit,i)
    y<-year.base+i
    if(!is.null(rs)){
      assign(".Random.seed", rs, envir = .GlobalEnv)
      rs <- NULL
    }
    sw<-getThisOrAve(fit$data$stockMeanWeight, y, ave.sw)
    cw<-getThisOrAve(fit$data$catchMeanWeight, y, ave.cw)
    mo<-getThisOrAve(fit$data$propMat, y, ave.mo)
    if(i>0){
      nmPreviousYear<-getThisOrAve(fit$data$natMor, y-1, ave.nm) *scaleM[i]
      nm<-getThisOrAve(fit$data$natMor, y, ave.nm) *scaleM[i+1]
    }else{
      nmPreviousYear<-getThisOrAve(fit$data$natMor, y-1, ave.nm)
      nm<-getThisOrAve(fit$data$natMor, y, ave.nm)
    }
    lf<-getThisOrAve(fit$data$landFrac, y, ave.lf)
    lw<-getThisOrAve(fit$data$landMeanWeight, y, ave.lw)
    pm<-getThisOrAve(fit$data$propM, y, ave.pm)
    pf<-getThisOrAve(fit$data$propF, y, ave.pf)
    if(useNMmodel){
      sim<-t(sapply(1:nrow(sim),
                    function(s){
                      yidx<-which(rnNM==as.integer(y))
                      thisnm<-matrix(exp(simLogNM[s,]),nrow=nrow(opar$logNM))[yidx,]
                      step(sim[s,], nm=thisnm, recpool=recpool, scale=1, inyear=(i==0))
                    }
      ))
    }else{
      sim <- t(apply(sim, 1, function(s)step(s, nm=nmPreviousYear, recpool=recpool, scale=1, inyear=(i==0))))
    }
    if(i!=0){
      cof <- coef(fit)
      if(deterministic){procVar<-procVar*0}
      if(!is.null(overwriteSelYears)){nn<-length(fit$conf$keyLogFsta[1,]); procVar[-c(1:nn),-c(1:nn)] <- 0}
      if(!is.null(customSel)){nn<-length(fit$conf$keyLogFsta[1,]); procVar[-c(1:nn),-c(1:nn)] <- 0}

      if(fit$conf$corFlag <3){
        sim <- sim + rmvnorm(nosim, mu=rep(0,nrow(procVar)), Sigma=procVar)
      }else if(fit$conf$corFlag ==3){
        k<-fit$conf$keyLogFsta[1,]
        k <- unique(k[k >= 0] + 1)
        nF <- length(k)
        simTmp <- rmvnorm(nosim, mu=rep(0,nrow(procVar)), Sigma=procVar)
        rhoF = 2/(1 + exp(-2*cof[names(cof)=="sepFlogitRho"])) -1
        sepFalpha = cof[names(cof)=="sepFalpha"]
        for(jj in 1:nosim){
          #Calculates previous U and V.
          V = mean(sim[jj,(dim(sim)[2]-nF+1):(dim(sim)[2])])
          U = sim[jj,(dim(sim)[2]-nF+1):(dim(sim)[2]-1)] -V

          #Extract new contributions to U and V
          Uepsilon = simTmp[jj,(dim(simTmp)[2]-nF+1):(dim(simTmp)[2]-1)]
          Vepsilon = simTmp[jj,dim(simTmp)[2]]

          #Calculates updated U and V
          Unew = rhoF[1]*U + Uepsilon+ sepFalpha[1:(length(sepFalpha)-1)]
          Vnew = rhoF[2]*V + Vepsilon+ sepFalpha[length(sepFalpha)]

          #Calculate F based on U and V
          sim[jj,(dim(simTmp)[2]-nF+1):(dim(simTmp)[2]-1)] = Unew + Vnew
          sim[jj,dim(simTmp)[2]] =  -sum(Unew) + Vnew #subtract sum(Unew) since U sum to 0 over all ages

          #Update process for N
          sim[jj,1:(dim(simTmp)[2]-nF)] = sim[jj,1:(dim(simTmp)[2]-nF)]   + simTmp[jj,1:(dim(simTmp)[2]-nF)] #Simulated N is not affected by the separable structure of F
        }
      }else{
        stop("forcast not implemented for the given corflag")
      }
    }


    if(!is.na(fscale[i+1])){
      sim<-t(apply(sim, 1, scaleF, scale=fscale[i+1]))
    }

    if(!is.na(fval[i+1])){
      curfbar<-median(apply(sim, 1, fbar))
      adj<-fval[i+1]/curfbar
      sim<-t(apply(sim, 1, scaleF, scale=adj))
    }

    if(!is.na(cwF[i+1])){
      if(missing(customWeights))stop("customWeights must be supplied when using the cwF option")
      curcwF<-median(apply(sim, 1, getCWF, w=customWeights))
      adj<-cwF[i+1]/curcwF
      sim<-t(apply(sim, 1, scaleF, scale=adj))
    }

    if(!is.na(catchval[i+1])){
      simtmp<-NA
      fun<-function(s){
        simtmp<<-t(apply(sim, 1, scaleF, scale=s))
        if(useCWmodel | useNMmodel){
          simcat<-sapply(1:nrow(simtmp), function(k){
            if(useCWmodel){
              thisy <- which(rnCW==y);
              thiscw <-matrix(exp(simLogCW[k,]),nrow=nrow(opar$logCW))[thisy,];
            }else{
              thiscw <- cw
            }
            if(useNMmodel){
              thisy <- which(rnNM==y);
              thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,];
            }else{
              thisnm <- nm
            }
            catch(simtmp[k,], nm=thisnm, cw=thiscw)}
          )
        }else{
          simcat<-apply(simtmp, 1, catch, nm=nm, cw=cw)
        }
        return(catchval[i+1]-median(simcat))
      }
      ff <- uniroot(fun, c(0,100))$root
      sim <- simtmp
    }

    if(!is.na(catchval.exact[i+1])){
      simtmp<-sim
      funk<-function(k){
        if(useCWmodel){
          thisy <- which(rnCW==y)
          thiscw <-matrix(exp(simLogCW[k,]),nrow=nrow(opar$logCW))[thisy,]
        }else{
          thiscw=cw
        }
        if(useNMmodel){
          thisy <- which(rnNM==y)
          thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,]
        }else{
          thisnm=nm
        }
        one<-function(s){
          simtmp[k,]<<-scaleF(sim[k,],s)
          simcat<-catch(simtmp[k,], nm=thisnm, cw=thiscw)
          return(catchval.exact[i+1]-simcat)
        }
        ff <- uniroot(one, c(0,100))$root
      }
      dd <- sapply(1:nrow(sim),funk)
      sim <- simtmp
    }

    if(!is.na(landval[i+1])){
      simtmp<-NA
      fun<-function(s){
        simtmp<<-t(apply(sim, 1, scaleF, scale=s))
        if(useNMmodel){
          simcat<-sapply(1:nrow(simtmp), function(k){
            thisy <- which(rnNM==y);
            thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,];
            catchFrac(simtmp[k,], nm=thisnm, w=lw, frac=lf)
          }
          )
        }else{
          simcat<-apply(simtmp, 1, catchFrac, nm=nm, w=lw, frac=lf)
        }
        return(landval[i+1]-median(simcat))
      }
      ff <- uniroot(fun, c(0,100))$root
      sim <- simtmp
    }

    if(!is.na(nextssb[i+1])){
      if((sum(pm)!=0) | (sum(pf)!=0))stop("nextssb option only available if SSB is calculated in the beginning of the year")
      simtmp<-NA
      rs<-.Random.seed
      fun<-function(s){
        assign(".Random.seed", rs, envir = .GlobalEnv)
        simtmp<<-t(apply(sim, 1, scaleF, scale=s))

        if(useNMmodel){
          simsim<-t(sapply(1:nrow(simtmp),
                           function(s){
                             yidx<-which(rnNM==as.integer(y))
                             thisnm<-matrix(exp(simLogNM[s,]),nrow=nrow(opar$logNM))[yidx,]
                             step(simtmp[s,], nm=thisnm, recpool=recpool, scale=1, inyear=(i==0))
                           }
          ))
        }else{
          simsim <- t(apply(simtmp, 1, function(s)step(s, nm=nm, recpool=recpool, scale=1, inyear=(i==0))))
        }
        if(i!=0){
          if(deterministic)procVar<-procVar*0
          simsim <- simsim + rmvnorm(nosim, mu=rep(0,nrow(procVar)), Sigma=procVar)
        }
        if(useSWmodel | useMOmodel | useNMmodel){
          ssbswmoi<-function(k){
            if(useSWmodel){
              thisy <- which(rnSW==y)
              thissw <-matrix(exp(simLogSw[k,]),nrow=nrow(opar$logSW))[thisy,]
            }else{
              thissw <-sw
            }
            if(useMOmodel){
              thisy <- which(rnMO==y)
              thismo <-matrix(plogis(simLogitMO[k,]),nrow=nrow(opar$logitMO))[thisy,]
            }else{
              thismo <-mo
            }
            if(useNMmodel){
              thisy <- which(rnNM==y)
              thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,]
            }else{
              thisnm <-nm
            }
            ssb(simsim[k,],nm=thisnm,sw=thissw,mo=thismo,pm=pm,pf=pf)
          }
          simnextssb <- sapply(1:nrow(simsim), ssbswmoi)
        }else{
          simnextssb <- apply(simsim, 1, ssb, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)
        }
        return(nextssb[i+1]-median(simnextssb))
      }
      ff <- uniroot(fun, c(0,100))$root
      sim <- simtmp
    }

    fbarsim <- apply(sim, 1, fbar)
    fbarLsim <- apply(sim, 1, fbarFrac, lf=lf)
    if(useCWmodel | useNMmodel){
      catchcwi<-function(k){
        if(useCWmodel){
          thisy <- which(rnCW==y)
          thiscw <-matrix(exp(simLogCW[k,]),nrow=nrow(opar$logCW))[thisy,]
        }else{
          thiscw <- cw
        }
        if(useNMmodel){
          thisy <- which(rnNM==y)
          thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,]
        }else{
          thisnm <- nm
        }
        catch(sim[k,], nm=thisnm, cw=thiscw)
      }
      catchsim <- sapply(1:nrow(sim), catchcwi)
    }else{
      catchsim <- apply(sim, 1, catch, nm=nm, cw=cw) #########################NB! Here we must use nm in simulation year##################
    }

    if(useNMmodel){
      landsim<-sapply(1:nrow(sim), function(k){
        thisy <- which(rnNM==y);
        thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,];
        catchFrac(sim[k,], nm=thisnm, w=lw, frac=lf)
      }
      )
    }else{
      landsim<-apply(sim, 1, catchFrac, nm=nm, w=lw, frac=lf)#########################NB! Here we must use nm in simulation year###################
    }
    catchatagesim <- apply(sim, 1, catchatage, nm=nm)######################NB! Here we must use nm in simulation year##############################
    if(useSWmodel | useMOmodel | useNMmodel){
      ssbswmoi<-function(k){
        if(useSWmodel){
          thisy <- which(rnSW==y)
          thissw <-matrix(exp(simLogSw[k,]),nrow=nrow(opar$logSW))[thisy,]
        }else{
          thissw <-sw
        }
        if(useMOmodel){
          thisy <- which(rnMO==y)
          thismo <-matrix(plogis(simLogitMO[k,]),nrow=nrow(opar$logitMO))[thisy,]
        }else{
          thismo <-mo
        }
        if(useNMmodel){
          thisy <- which(rnNM==y)
          thisnm <-matrix(exp(simLogNM[k,]),nrow=nrow(opar$logNM))[thisy,]
        }else{
          thisnm <-nm
        }
        ssb(sim[k,],nm=thisnm,sw=thissw,mo=thismo,pm=pm,pf=pf)
      }
      ssbsim <- sapply(1:nrow(sim), ssbswmoi)
    }else{
      ssbsim <- apply(sim, 1, ssb, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)######################NB! Here we must use nm in simulation year##############################
    }
    tsbsim <- apply(sim, 1, tsb, sw=sw)
    if(lagR){
      recsim <- exp(sim[,2])
    }else{
      recsim <- exp(sim[,1])
    }
    cwFsim <- rep(NA,nrow(sim))
    if(!missing(customWeights)){
      cwFsim <- apply(sim, 1, getCWF, w=customWeights)
    }
    simlist[[i+1]] <- list(sim=sim, fbar=fbarsim, catch=catchsim, ssb=ssbsim, rec=recsim,
                           cwF=cwFsim, catchatage=catchatagesim, land=landsim, fbarL=fbarLsim, tsb=tsbsim, year=y)
  }

  attr(simlist, "fit")<-fit

  collect <- function(x){
    quan <- quantile(x, c(.50,.025,.975))
    c(median=quan[1], low=quan[2], high=quan[3])
  }

  fbar <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbar))),3)
  fbarL <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbarL))),3)
  rec <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$rec))))
  ssb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssb))))
  tsb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$tsb))))
  catch <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$catch))))
  land <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$land))))
  caytable<-round(do.call(rbind, lapply(simlist, function(xx)apply(xx$catchatage,1,collect))))
  tab <- cbind(fbar,rec,ssb,catch)
  if(splitLD){
    tab<-cbind(tab,fbarL,fbar-fbarL,land,catch-land)
  }
  if(addTSB){
    tab<-cbind(tab,tsb)
  }
  if(!missing(customWeights)) tab <- cbind(tab,cwF=round(do.call(rbind, lapply(simlist, function(xx)collect(xx$cwF))),3))
  rownames(tab) <- unlist(lapply(simlist, function(xx)xx$year))
  nam <- c("median","low","high")
  basename<-c("fbar:","rec:","ssb:","catch:")
  if(splitLD){
    basename<-c(basename,"fbarL:","fbarD:","Land:","Discard:")
  }
  if(addTSB){
    basename<-c(basename,"tsb:")
  }
  if(!missing(customWeights)){
    basename<-c(basename,"cwF:")
  }
  colnames(tab)<-paste0(rep(basename, each=length(nam)), nam)
  attr(simlist, "tab")<-tab
  shorttab<-t(tab[,grep("median",colnames(tab))])
  rownames(shorttab)<-sub(":median","",paste0(label,if(!is.null(label))":",rownames(shorttab)))
  attr(simlist, "shorttab")<-shorttab
  attr(simlist, "label") <- label
  attr(simlist, "caytable")<-caytable
  class(simlist) <- "samforecast"
  if(!savesim){
    simlistsmall<-lapply(simlist, function(x)list(year=x$year))
    attributes(simlistsmall)<-attributes(simlist)
    return(simlistsmall)
  }else{
    return(simlist)
  }
}





quantileRetro = function(rets){

  fit = attributes(rets[[1]])$fit
  lagF = max(fit$data$years) -max(as.numeric(rownames(getFleet(fit,1))))

  if(lagF>1){
    stop("Catch not given in terminal year or terminal year-1")
  }

  lag = c(0,0,lagF)

  tableOfInterest <- function(fit) {
    ret <- cbind(rectable(fit)[, 1], ssbtable(fit)[,1], fbartable(fit)[,1])
    colnames(ret) <- c(paste("R(age ", fit$conf$minAge, ")", sep = ""), "SSB",
                       paste("Fbar(", fit$conf$fbarRange[1], "-",
                             fit$conf$fbarRange[2], ")", sep = ""))
    ret
  }

  biasMat = list()
  nsim = length(rets)
  nya = length(rets[[1]])
  for(i in 1:nsim){
    fits = rets[[i]]
    ref <- tableOfInterest(attr(fits, "fit"))
    ret <- lapply(fits, tableOfInterest)

    bias <- lapply(ret, function(x) {
      y <- rownames(x)[nrow(x) - lag]
      kk = rep(0,3)
      for(j in 1:3){
        kk[j] = (x[rownames(x) == y[j],j ]-ref[rownames(ref) == y[j],j] )/ref[rownames(ref) == y[j],j]
      }
      kk
    })
    biasMat[[i]] = do.call(rbind, bias)
  }

  lQ = matrix(0,nrow = nya,ncol = 3)
  uQ = matrix(0,nrow = nya,ncol = 3)

  biasArray = array(as.numeric(unlist(biasMat)), dim=c(nya, 3, nsim))


  for(y in 1:nya){
    q <- apply(biasArray[y,,], 1, quantile, probs = c(0.025,0.975),  na.rm = TRUE)
    lQ[y,] = q[1,]
    uQ[y,] = q[2,]
  }

  lQ = as.data.frame(lQ); uQ = as.data.frame(uQ);
  colnames(lQ) = c("R","SSB","F")
  colnames(uQ) = c("R","SSB","F")

  return(list(lQ = lQ, uQ = uQ))
}

