#Construct confidence intervals of Mohn's rho. Note that NCCcod.R and NEAcod.R needs to be ran first to save Mohn's rho samples

library(stockassessment)
##################------------------------------------------NEA cod----------------------------------##################################
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

ar1 = FALSE
nya = 5
if(!ar1){
  ret = retro(fit,year = nya)
  png("figures/mohnNEAcod.png")
  mohT = read.table("results/NEACodOfficial.txt", header = TRUE)
}else{
  conf$corFlag=2
  par = defpar(dat,conf)
  fitAR1 = sam.fit(dat,conf,par)
  ret = retro(fitAR1,year = nya)
  png("figures/mohnNEAcodAR1.png")
  mohT = read.table("results/NEACodAR1.txt", header = TRUE)
}
p = 0.05
nya = 5
lagF = max(fit$data$years) -max(as.numeric(rownames(getFleet(fit,1))))
qR = quantile(mohT[,1], c(p/2,1-p/2))
qSSB = quantile(mohT[,2], c(p/2,1-p/2))
qF = quantile(mohT[,3], c(p/2,1-p/2))

#ylim = c(min(mohT), max(mohT))
ylim = c(-0.4, 0.4)
library(plotrix)
plotCI(3,y = median(mohT[,1]), li = qR[1], ui = qR[2], ylim = ylim, xlim = c(0.5,3.3),
       main = "",
       xlab = "",
       ylab = "Mohn's rho",
       xaxt='n', lwd = 2,cex.axis = 1.5,cex.lab = 1.5)
for(i in 1:dim(mohT)[1]){
#  if(mohT[i,1] > qR[2] | mohT[i,1] < qR[1]) points(3,mohT[i,1], cex = 0.5)
}

plotCI(1,y = median(mohT[,2]), li = qSSB[1], ui = qSSB[2], ylim = ylim, lwd = 2, add = TRUE)
for(i in 2:dim(mohT)[1]){
#  if(mohT[i,2] > qSSB[2] | mohT[i,2] < qSSB[1]) points(1,mohT[i,2], cex = 0.5)
}

plotCI(2,y = median(mohT[,3]), li = qF[1], ui = qF[2], ylim = ylim, lwd = 2, add = TRUE)
for(i in 2:dim(mohT)[1]){
#  if(mohT[i,3] > qF[2] | mohT[i,3] < qF[1]) points(2,mohT[i,3], cex = 0.5)
}

abline(h = 0)


axis(1, at=1:3, labels=c("SSB","F-bar","R"), las=1,cex.axis = 1.5)
mohRObs = mohn(ret)[1]
mohSSBObs = mohn(ret)[2]
mohFObs = mohn(ret,lag = lagF)[3]
mohObs = c(mohRObs,mohSSBObs,mohFObs)
points(mohObs[c(2,3,1)], col = 'red', cex = 1, pch = 16)


pp = rep(0,3)
for(i in 1:3){
  pp[i] = round(min(mean(mohT[,i]>mohObs[i]),mean(mohT[,i]<mohObs[i]))*2,2)
}
if(pp[2]<p) col = 'black' else col = "black"
text(0.7, 0.02, paste0("P:  ",pp[2]),
     cex=1.3, pos=3,col=col)
if(pp[3]<p) col = 'black' else col = "black"
text(1.7, 0.02, paste0("P:  ",pp[3]),
     cex=1.3, pos=3,col=col)
if(pp[1]<p) col = 'black' else col = "black"
text(2.7, 0.02, paste0("P:  ",pp[1]),
     cex=1.3, pos=3,col=col)


points(c(0.6,1.4),c(0.2,0.2), type = 'l', lty = 2,col = 'green',lwd = 2)
points(c(0.6,1.4),c(-0.15,-0.15), type = 'l', lty = 2,col = 'green',lwd = 2)


title("Mohn's rho distribution NEA cod", cex.main = 2)

if(ar1){
  ret2 = retro(fit,year = nya)
  mohT2 = read.table("results/NEACodOfficial.txt", header = TRUE)

  qR = quantile(mohT2[,1], c(p/2,1-p/2))
  qSSB = quantile(mohT2[,2], c(p/2,1-p/2))
  qF = quantile(mohT2[,3], c(p/2,1-p/2))

  ylim = c(-0.4, 0.4)
  plotCI(3.1,y = median(mohT2[,1]), li = qR[1], ui = qR[2], add = TRUE,col = 'grey',lwd = 2)
  for(i in 1:dim(mohT2)[1]){
#    if(mohT2[i,1] > qR[2] | mohT2[i,1] < qR[1]) points(3.1,mohT2[i,1], cex = 0.5,col = 'grey')
  }

  plotCI(1.1,y = median(mohT2[,2]), li = qSSB[1], ui = qSSB[2], ylim = ylim, lwd = 2, add = TRUE,col = 'grey')
  for(i in 2:dim(mohT2)[1]){
#    if(mohT2[i,2] > qSSB[2] | mohT2[i,2] < qSSB[1]) points(1.1,mohT2[i,2], cex = 0.5,col = 'grey')
  }

  plotCI(2.1,y = median(mohT2[,3]), li = qF[1], ui = qF[2], ylim = ylim, lwd = 2, add = TRUE,col = 'grey')
  for(i in 2:dim(mohT2)[1]){
#    if(mohT2[i,3] > qF[2] | mohT2[i,3] < qF[1]) points(2.1,mohT2[i,3], cex = 0.5,col = 'grey')
  }


  mohRObs = mohn(ret2)[1]
  mohSSBObs = mohn(ret2)[2]
  mohFObs = mohn(ret2,lag = lagF)[3]
  mohObs = c(mohRObs,mohSSBObs,mohFObs)
  points(c(1.1,2.1,3.1),mohObs[c(2,3,1)], col = 'grey', cex = 1, pch = 16)

}

dev.off()








##################------------------------------------------coastal cod----------------------------------##################################


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


conf = loadConf(dat,file = "data/coastalCod/conf.cfg")
par = defpar(dat,conf)
fit = sam.fit(dat,conf,par)
nya = 5
ret = retro(fit,year = nya)
mohT = read.table("results/costalCodOfficial.txt", header = TRUE)

p = 0.05
lagF = max(fit$data$years) -max(as.numeric(rownames(getFleet(fit,1))))


png("figures/mohnCostalcod.png")
qR = quantile(mohT[,1], c(p/2,1-p/2))
qSSB = quantile(mohT[,2], c(p/2,1-p/2))
qF = quantile(mohT[,3], c(p/2,1-p/2))

#ylim = c(min(mohT), max(mohT))
ylim = c(-0.4, 0.4)
library(plotrix)
plotCI(3,y = median(mohT[,1]), li = qR[1], ui = qR[2], ylim = ylim, xlim = c(0.5,3.3),
       main = "",
       xlab = "",
       ylab = "Mohn's rho",
       xaxt='n', lwd = 2,cex.axis = 1.5,cex.lab = 1.5)
for(i in 1:dim(mohT)[1]){
#  if(mohT[i,1] > qR[2] | mohT[i,1] < qR[1]) points(3,mohT[i,1], cex = 0.5)
}

plotCI(1,y = median(mohT[,2]), li = qSSB[1], ui = qSSB[2], ylim = ylim, lwd = 2, add = TRUE)
for(i in 2:dim(mohT)[1]){
#  if(mohT[i,2] > qSSB[2] | mohT[i,2] < qSSB[1]) points(1,mohT[i,2], cex = 0.5)
}

plotCI(2,y = median(mohT[,3]), li = qF[1], ui = qF[2], ylim = ylim, lwd = 2, add = TRUE)
for(i in 2:dim(mohT)[1]){
#  if(mohT[i,3] > qF[2] | mohT[i,3] < qF[1]) points(2,mohT[i,3], cex = 0.5)
}

abline(h = 0)

axis(1, at=1:3, labels=c("SSB","F-bar","R"), las=1,cex.axis = 1.5)
mohRObs = mohn(ret)[1]
mohSSBObs = mohn(ret)[2]
mohFObs = mohn(ret,lag = lagF)[3]
mohObs = c(mohRObs,mohSSBObs,mohFObs)
points(mohObs[c(2,3,1)], col = 'red', cex = 1, pch = 16)


pp = rep(0,3)
for(i in 1:3){
  pp[i] = round(min(mean(mohT[,i]>mohObs[i]),mean(mohT[,i]<mohObs[i]))*2,2)
}
if(pp[2]<p) col = 'black' else col = "black"
text(0.7, 0.02, paste0("P:  ",pp[2]),
     cex=1.3, pos=3,col=col)
if(pp[3]<p) col = 'black' else col = "black"
text(1.7, 0.02, paste0("P:  ",pp[3]),
     cex=1.3, pos=3,col=col)
if(pp[1]<p) col = 'black' else col = "black"
text(2.7, 0.02, paste0("P:  ",pp[1]),
     cex=1.3, pos=3,col=col)


points(c(0.6,1.4),c(0.2,0.2), type = 'l', lty = 2,col = 'green',lwd = 2)
points(c(0.6,1.4),c(-0.15,-0.15), type = 'l', lty = 2,col = 'green',lwd = 2)

title("Mohn's rho distribution NCC cod", cex.main = 2)
dev.off()



