#Let "rets" be the list with simulated retros and "ret" be the retro
fit = attributes(rets[[1]])$fit
lagF = max(fit$data$years) -max(as.numeric(rownames(getFleet(fit,1))))
xlim = c(2010,2020)

fitFirst = ret[[length(ret)]]

png(filename = "figures/FsimAR1.png",width = 1000,height = 1200)
par(mfrow = c(5,5),mar = c(1, 3, 3, 1))

ff = attributes(ret)$fit
ylim = c(0, 1.2)
fbarplot(ret,main = paste0("Fbar obtained"),xlim = xlim,ylim = ylim)
for(ii in 1:24){
  fbar = calcFbar(attributes(rets)$simF[[ii]], ageRange = attributes(rets)$fit$conf$fbarRange)
  tmp = (dim(attributes(rets)$fit$pl$logN)[2]-nya):dim(attributes(rets)$fit$pl$logN)[2]
  ff = attributes(rets[[ii]])$fit
#  ylim = c(0, max(fbartable(ff)[(dim(fbartable(ff))[1]-10):(dim(fbartable(ff))[1]-1),3])+0.1)

  fbarplot(rets[[ii]],main = paste0("Fbar, simulation ",ii),xlim = xlim,ylim = ylim)
  points(max(fitFirst$data$years) + 0:(length(tmp)-1-lagF),fbar[1:(nya +1-lagF)],type = 'l', col = 2)
}

dev.off()


#SSB
par(mfrow = c(5,5),mar = c(1, 3, 3, 1))

ff = attributes(ret)$fit
ylim = c(0, max(ssbtable(ff)[(dim(ssbtable(ff))[1]-10):dim(ssbtable(ff))[1],3]))
ssbplot(ret,main = paste0("SSB, obtained ",1),xlim = xlim,ylim = ylim)

for(ii in 1:24){
  ff = attributes(rets[[ii]])$fit
  ylim = c(0, max(ssbtable(ff)[(dim(ssbtable(ff))[1]-10):dim(ssbtable(ff))[1],3]))
  tmp = (dim(attributes(rets)$fit$pl$logN)[2]-nya):dim(attributes(rets)$fit$pl$logN)[2]
  ssb = calcSSB(attributes(rets)$simN[[ii]], attributes(rets)$fit$data$stockMeanWeight[tmp,], attributes(rets)$fit$data$propMat[tmp,],
                attributes(rets)$fit$data$natMor[tmp,], attributes(rets)$fit$data$propF[tmp,],
                attributes(rets)$fit$data$propM[tmp,], attributes(rets)$simF[[ii]])

  ssbplot(rets[[ii]],main = paste0("SSB, simulation ",ii),xlim = xlim,ylim = ylim)
  points(max(fitFirst$data$years) + 0:(length(tmp)-1),ssb,type = 'l', col = 2)
}


