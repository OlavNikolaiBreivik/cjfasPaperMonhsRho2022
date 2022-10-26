#Construct strength plots. Note that NCCcod.R and NEAcod.R needs to be ran first to save Mohn's rho samples

##################------------------------------------------NEA cod----------------------------------##################################
#R
png(filename = "figures/neaCodStrengthR.png",width = 700)
mohT = read.table("results/NEACodOfficial.txt", header = TRUE)
mohM3 = read.table("results/NEACodM3.txt", header = TRUE)
mohM05 = read.table("results/NEACodM05.txt", header = TRUE)
mohQ3 = read.table("results/NEACodQ3.txt", header = TRUE)
mohC05 = read.table("results/NEACodC05.txt", header = TRUE)
mohSdF3 = read.table("results/NEACodSdF3.txt", header = TRUE)
mohSdObsI3 = read.table("results/NEACodSdObsI3.txt", header = TRUE)
mohSdObsC3 = read.table("results/NEACodSdObsC3.txt", header = TRUE)
mohCorObs = read.table("results/NEACodCorObs.txt", header = TRUE)

moh = data.frame(mohT[,1])
moh$M3 = mohM3[,1]
moh$M05 = mohM05[,1]
moh$Q3 = mohQ3[,1]
moh$C05 = mohC05[,1]
moh$SdF3 = mohSdF3[,1]
moh$SdObsI3 = mohSdObsI3[,1]
moh$SdObsC3 = mohSdObsC3[,1]
moh$CorObs = mohCorObs[,1]


library(plotrix)
p = 0.05
i = 1
par(mar=c(4,6,4,4))
plotCI(i,median(moh[,2]), li = quantile(moh[,2],p/2), ui =quantile(moh[,2],1-p/2), xlim = c(0.3,8), ylim = c(-0.5,0.8),
       ylab = "Mohn's rho", xlab = "Scenario", main = "Test power, recruitment NEA cod", xaxt= "n", cex.axis = 1.5,cex.lab = 2.2,
       cex.main = 3, lwd = 2)
i = 2
plotCI(i,median(moh[,3]), li = quantile(moh[,3],p/2), ui =quantile(moh[,3],1-p/2), add = TRUE, lwd = 2)
i = 3
plotCI(i,median(moh[,4]), li = quantile(moh[,4],p/2), ui =quantile(moh[,4],1-p/2), add = TRUE, lwd = 2)
i = 4
plotCI(i,median(moh[,5]), li = quantile(moh[,5],p/2), ui =quantile(moh[,5],1-p/2), add = TRUE, lwd = 2)
i = 5
plotCI(i,median(moh[,6]), li = quantile(moh[,6],p/2), ui =quantile(moh[,6],1-p/2), add = TRUE, lwd = 2)
i = 6
plotCI(i,median(moh[,7]), li = quantile(moh[,7],p/2), ui =quantile(moh[,7],1-p/2), add = TRUE, lwd = 2)
i = 7
plotCI(i,median(moh[,8]), li = quantile(moh[,8],p/2), ui =quantile(moh[,8],1-p/2), add = TRUE, lwd = 2)
i = 8
plotCI(i,median(moh[,9]), li = quantile(moh[,9],p/2), ui =quantile(moh[,9],1-p/2), add = TRUE, lwd = 2)


abline(h = 0)
abline(h = quantile(moh[,1],p/2), col = 'green' ,lty = 2, lwd = 2)
abline(h = quantile(moh[,1],1-p/2), col = 'green' ,lty = 2, lwd = 2)

S = round(mean(moh[,1]<quantile(moh[,1],p/2) | moh[,1]>quantile(moh[,1],1-p/2)),2)
SRM3 = round(mean(moh[,2]<quantile(moh[,1],p/2) | moh[,2]>quantile(moh[,1],1-p/2)),2)
SRM05 = round(mean(moh[,3]<quantile(moh[,1],p/2) | moh[,3]>quantile(moh[,1],1-p/2)),2)
SRQ3 = round(mean(moh[,4]<quantile(moh[,1],p/2) | moh[,4]>quantile(moh[,1],1-p/2)),2)
SRC05 = round(mean(moh[,5]<quantile(moh[,1],p/2) | moh[,5]>quantile(moh[,1],1-p/2)),2)
SRSdF3 = round(mean(moh[,6]<quantile(moh[,1],p/2) | moh[,6]>quantile(moh[,1],1-p/2)),2)
SRSdObsI3 = round(mean(moh[,7]<quantile(moh[,1],p/2) | moh[,7]>quantile(moh[,1],1-p/2)),2)
SRSdObsC3 = round(mean(moh[,8]<quantile(moh[,1],p/2) | moh[,8]>quantile(moh[,1],1-p/2)),2)
SRCorObs = round(mean(moh[,9]<quantile(moh[,1],p/2) | moh[,9]>quantile(moh[,1],1-p/2)),2)

text(0.5,median(moh[,2])-0.4, paste0("S:",SRM3),
     cex=1.3, pos=3)
text(1.5,median(moh[,3]), paste0("S:",SRM05),
     cex=1.3, pos=3)
text(2.5,median(moh[,4]), paste0("S:",SRQ3),
     cex=1.3, pos=3)
text(3.5,median(moh[,5]), paste0("S:",SRC05),
     cex=1.3, pos=3)
text(4.5,median(moh[,6]), paste0("S:",SRSdF3),
     cex=1.3, pos=3)
text(5.5,median(moh[,7]), paste0("S:",SRSdObsI3),
     cex=1.3, pos=3)
text(6.5,median(moh[,8]), paste0("S:",SRSdObsC3),
     cex=1.3, pos=3)
text(7.5,median(moh[,9]), paste0("S:",SRCorObs),
     cex=1.3, pos=3)

axis(side = 1, at = 1:8, labels = c("M3","M/3","Q3","C/3","sdF", "sdI", "sdC", "corObs"), tck = -0.01, cex.axis = 1.5)
dev.off()

#SSB
png(filename = "figures/neaCodStrengthSSB.png",width = 700)

moh = data.frame(mohT[,2])
moh$M3 = mohM3[,2]
moh$M05 = mohM05[,2]
moh$Q3 = mohQ3[,2]
moh$C05 = mohC05[,2]
moh$SdF3 = mohSdF3[,2]
moh$SdObsI3 = mohSdObsI3[,2]
moh$SdObsC3 = mohSdObsC3[,2]
moh$CorObs = mohCorObs[,2]

library(plotrix)
p = 0.05
i = 1
par(mar=c(4,6,4,4))
plotCI(i,median(moh[,2]), li = quantile(moh[,2],p/2), ui =quantile(moh[,2],1-p/2), xlim = c(0.3,8), ylim = c(-0.5,0.8),
       ylab = "Mohn's rho", xlab = "Scenario", main = "Test power, SSB NEA cod", xaxt= "n", cex.axis = 1.5,cex.lab = 2.2,
       cex.main = 3, lwd = 2)
i = 2
plotCI(i,median(moh[,3]), li = quantile(moh[,3],p/2), ui =quantile(moh[,3],1-p/2), add = TRUE, lwd = 2)
i = 3
plotCI(i,median(moh[,4]), li = quantile(moh[,4],p/2), ui =quantile(moh[,4],1-p/2), add = TRUE, lwd = 2)
i = 4
plotCI(i,median(moh[,5]), li = quantile(moh[,5],p/2), ui =quantile(moh[,5],1-p/2), add = TRUE, lwd = 2)
i = 5
plotCI(i,median(moh[,6]), li = quantile(moh[,6],p/2), ui =quantile(moh[,6],1-p/2), add = TRUE, lwd = 2)
i = 6
plotCI(i,median(moh[,7]), li = quantile(moh[,7],p/2), ui =quantile(moh[,7],1-p/2), add = TRUE, lwd = 2)
i = 7
plotCI(i,median(moh[,8]), li = quantile(moh[,8],p/2), ui =quantile(moh[,8],1-p/2), add = TRUE, lwd = 2)
i = 8
plotCI(i,median(moh[,9]), li = quantile(moh[,9],p/2), ui =quantile(moh[,9],1-p/2), add = TRUE, lwd = 2)
abline(h = 0)
abline(h = quantile(moh[,1],p/2), col = 'green' ,lty = 2)
abline(h = quantile(moh[,1],1-p/2), col = 'green' ,lty = 2)

S = round(mean(moh[,1]<quantile(moh[,1],p/2) | moh[,1]>quantile(moh[,1],1-p/2)),2)
SRM3 = round(mean(moh[,2]<quantile(moh[,1],p/2) | moh[,2]>quantile(moh[,1],1-p/2)),2)
SRM05 = round(mean(moh[,3]<quantile(moh[,1],p/2) | moh[,3]>quantile(moh[,1],1-p/2)),2)
SRQ3 = round(mean(moh[,4]<quantile(moh[,1],p/2) | moh[,4]>quantile(moh[,1],1-p/2)),2)
SRC05 = round(mean(moh[,5]<quantile(moh[,1],p/2) | moh[,5]>quantile(moh[,1],1-p/2)),2)
SRSdF3 = round(mean(moh[,6]<quantile(moh[,1],p/2) | moh[,6]>quantile(moh[,1],1-p/2)),2)
SRSdObsI3 = round(mean(moh[,7]<quantile(moh[,1],p/2) | moh[,7]>quantile(moh[,1],1-p/2)),2)
SRSdObsC3 = round(mean(moh[,8]<quantile(moh[,1],p/2) | moh[,8]>quantile(moh[,1],1-p/2)),2)
SRCorObs = round(mean(moh[,9]<quantile(moh[,1],p/2) | moh[,9]>quantile(moh[,1],1-p/2)),2)

text(0.5,median(moh[,2]), paste0("S:",SRM3),
     cex=1.3, pos=3)
text(1.5,median(moh[,3]), paste0("S:",SRM05),
     cex=1.3, pos=3)
text(2.5,median(moh[,4]), paste0("S:",SRQ3),
     cex=1.3, pos=3)
text(3.5,median(moh[,5]), paste0("S:",SRC05),
     cex=1.3, pos=3)
text(4.5,median(moh[,6]), paste0("S:",SRSdF3),
     cex=1.3, pos=3)
text(5.5,median(moh[,7]), paste0("S:",SRSdObsI3),
     cex=1.3, pos=3)
text(6.5,median(moh[,8]), paste0("S:",SRSdObsC3),
     cex=1.3, pos=3)
text(7.5,median(moh[,9]), paste0("S:",SRCorObs),
     cex=1.3, pos=3)

axis(side = 1, at = 1:8, labels = c("M3","M/3","Q3","C/3","sdF", "sdI", "sdC", "corObs"), tck = -0.01, cex.axis = 1.5)
dev.off()

#F-bar
png(filename = "figures/neaCodStrengthF.png",width = 700)

moh = data.frame(mohT[,3])
moh$M3 = mohM3[,3]
moh$M05 = mohM05[,3]
moh$Q3 = mohQ3[,3]
moh$C05 = mohC05[,3]
moh$SdF3 = mohSdF3[,3]
moh$SdObsI3 = mohSdObsI3[,3]
moh$SdObsC3 = mohSdObsC3[,3]
moh$CorObs = mohCorObs[,3]
library(plotrix)
p = 0.05
par(mar=c(4,6,4,4))
i = 1
plotCI(i,median(moh[,2]), li = quantile(moh[,2],p/2), ui =quantile(moh[,2],1-p/2), xlim = c(0.3,8), ylim = c(-0.5,0.8),
       ylab = "Mohn's rho", xlab = "Scenario", main = "Test power, F-bar NEA cod", xaxt= "n", cex.axis = 1.5,cex.lab = 2.2,
       cex.main = 3, lwd = 2)
i = 2
plotCI(i,median(moh[,3]), li = quantile(moh[,3],p/2), ui =quantile(moh[,3],1-p/2), add = TRUE, lwd = 2)
i = 3
plotCI(i,median(moh[,4]), li = quantile(moh[,4],p/2), ui =quantile(moh[,4],1-p/2), add = TRUE, lwd = 2)
i = 4
plotCI(i,median(moh[,5]), li = quantile(moh[,5],p/2), ui =quantile(moh[,5],1-p/2), add = TRUE, lwd = 2)
i = 5
plotCI(i,median(moh[,6]), li = quantile(moh[,6],p/2), ui =quantile(moh[,6],1-p/2), add = TRUE, lwd = 2)
i = 6
plotCI(i,median(moh[,7]), li = quantile(moh[,7],p/2), ui =quantile(moh[,7],1-p/2), add = TRUE, lwd = 2)
i = 7
plotCI(i,median(moh[,8]), li = quantile(moh[,8],p/2), ui =quantile(moh[,8],1-p/2), add = TRUE, lwd = 2)
i = 8
plotCI(i,median(moh[,9]), li = quantile(moh[,9],p/2), ui =quantile(moh[,9],1-p/2), add = TRUE, lwd = 2)
abline(h = 0)
abline(h = quantile(moh[,1],p/2), col = 'green' ,lty = 2)
abline(h = quantile(moh[,1],1-p/2), col = 'green' ,lty = 2)

S = round(mean(moh[,1]<quantile(moh[,1],p/2) | moh[,1]>quantile(moh[,1],1-p/2)),2)
SRM3 = round(mean(moh[,2]<quantile(moh[,1],p/2) | moh[,2]>quantile(moh[,1],1-p/2)),2)
SRM05 = round(mean(moh[,3]<quantile(moh[,1],p/2) | moh[,3]>quantile(moh[,1],1-p/2)),2)
SRQ3 = round(mean(moh[,4]<quantile(moh[,1],p/2) | moh[,4]>quantile(moh[,1],1-p/2)),2)
SRC05 = round(mean(moh[,5]<quantile(moh[,1],p/2) | moh[,5]>quantile(moh[,1],1-p/2)),2)
SRSdF3 = round(mean(moh[,6]<quantile(moh[,1],p/2) | moh[,6]>quantile(moh[,1],1-p/2)),2)
SRSdObsI3 = round(mean(moh[,7]<quantile(moh[,1],p/2) | moh[,7]>quantile(moh[,1],1-p/2)),2)
SRSdObsC3 = round(mean(moh[,8]<quantile(moh[,1],p/2) | moh[,8]>quantile(moh[,1],1-p/2)),2)
SRCorObs = round(mean(moh[,9]<quantile(moh[,1],p/2) | moh[,9]>quantile(moh[,1],1-p/2)),2)

text(0.5,median(moh[,2]), paste0("S:",SRM3),
     cex=1.3, pos=3)
text(1.5,median(moh[,3]), paste0("S:",SRM05),
     cex=1.3, pos=3)
text(2.5,median(moh[,4]), paste0("S:",SRQ3),
     cex=1.3, pos=3)
text(3.5,median(moh[,5]), paste0("S:",SRC05),
     cex=1.3, pos=3)
text(4.5,median(moh[,6]), paste0("S:",SRSdF3),
     cex=1.3, pos=3)
text(5.5,median(moh[,7]), paste0("S:",SRSdObsI3),
     cex=1.3, pos=3)
text(6.5,median(moh[,8]), paste0("S:",SRSdObsC3),
     cex=1.3, pos=3)
text(7.5,median(moh[,9]), paste0("S:",SRCorObs),
     cex=1.3, pos=3)

axis(side = 1, at = 1:8, labels = c("M3","M/3","Q3","C/3","sdF", "sdI", "sdC", "corObs"), tck = -0.01, cex.axis = 1.5)

dev.off()








##################------------------------------------------Costal cod----------------------------------##################################
##################------------------------------------------Costal cod----------------------------------##################################
#R
png(filename = "figures/costalCodStrengthR.png",width = 700)
mohT = read.table("results/costalCodOfficial.txt", header = TRUE)
mohM3 = read.table("results/costalCodM3.txt", header = TRUE)
mohM05 = read.table("results/costalCodM05.txt", header = TRUE)
mohQ3 = read.table("results/costalCodQ3.txt", header = TRUE)
mohC05 = read.table("results/costalCodC05.txt", header = TRUE)
mohSdF3 = read.table("results/costalCodSdF3.txt", header = TRUE)
mohSdObsI3 = read.table("results/costalCodSdObsI3.txt", header = TRUE)
mohSdObsC3 = read.table("results/costalCodSdObsC3.txt", header = TRUE)
mohCorObs = read.table("results/costalCodCorObs.txt", header = TRUE)

moh = data.frame(mohT[,1])
moh$M3 = mohM3[,1]
moh$M05 = mohM05[,1]
moh$Q3 = mohQ3[,1]
moh$C05 = mohC05[,1]
moh$SdF3 = mohSdF3[,1]
moh$SdObsI3 = mohSdObsI3[,1]
moh$SdObsC3 = mohSdObsC3[,1]
moh$CorObs = mohCorObs[,1]


library(plotrix)
p = 0.05
par(mar=c(4,6,4,4))
i = 1
plotCI(i,median(moh[,2]), li = quantile(moh[,2],p/2), ui =quantile(moh[,2],1-p/2), xlim = c(0.3,8), ylim = c(-0.7,1),
       ylab = "Mohn's rho", xlab = "Scenario", main = "Test power, recruitment NCC cod", xaxt= "n", cex.axis = 1.5,cex.lab = 2.2,
       cex.main = 3, lwd = 2)
i = 2
plotCI(i,median(moh[,3]), li = quantile(moh[,3],p/2), ui =quantile(moh[,3],1-p/2), add = TRUE, lwd = 2)
i = 3
plotCI(i,median(moh[,4]), li = quantile(moh[,4],p/2), ui =quantile(moh[,4],1-p/2), add = TRUE, lwd = 2)
i = 4
plotCI(i,median(moh[,5]), li = quantile(moh[,5],p/2), ui =quantile(moh[,5],1-p/2), add = TRUE, lwd = 2)
i = 5
plotCI(i,median(moh[,6]), li = quantile(moh[,6],p/2), ui =quantile(moh[,6],1-p/2), add = TRUE, lwd = 2)
i = 6
plotCI(i,median(moh[,7]), li = quantile(moh[,7],p/2), ui =quantile(moh[,7],1-p/2), add = TRUE, lwd = 2)
i = 7
plotCI(i,median(moh[,8]), li = quantile(moh[,8],p/2), ui =quantile(moh[,8],1-p/2), add = TRUE, lwd = 2)
i = 8
plotCI(i,median(moh[,9]), li = quantile(moh[,9],p/2), ui =quantile(moh[,9],1-p/2), add = TRUE, lwd = 2)


abline(h = 0)
abline(h = quantile(moh[,1],p/2), col = 'green' ,lty = 2)
abline(h = quantile(moh[,1],1-p/2), col = 'green' ,lty = 2)

S = round(mean(moh[,1]<quantile(moh[,1],p/2) | moh[,1]>quantile(moh[,1],1-p/2)),2)
SRM3 = round(mean(moh[,2]<quantile(moh[,1],p/2) | moh[,2]>quantile(moh[,1],1-p/2)),2)
SRM05 = round(mean(moh[,3]<quantile(moh[,1],p/2) | moh[,3]>quantile(moh[,1],1-p/2)),2)
SRQ3 = round(mean(moh[,4]<quantile(moh[,1],p/2) | moh[,4]>quantile(moh[,1],1-p/2)),2)
SRC05 = round(mean(moh[,5]<quantile(moh[,1],p/2) | moh[,5]>quantile(moh[,1],1-p/2)),2)
SRSdF3 = round(mean(moh[,6]<quantile(moh[,1],p/2) | moh[,6]>quantile(moh[,1],1-p/2)),2)
SRSdObsI3 = round(mean(moh[,7]<quantile(moh[,1],p/2) | moh[,7]>quantile(moh[,1],1-p/2)),2)
SRSdObsC3 = round(mean(moh[,8]<quantile(moh[,1],p/2) | moh[,8]>quantile(moh[,1],1-p/2)),2)
SRCorObs = round(mean(moh[,9]<quantile(moh[,1],p/2) | moh[,9]>quantile(moh[,1],1-p/2)),2)

text(0.5,min(0.7,median(moh[,2])), paste0("S:",SRM3),
     cex=1.3, pos=3)
text(1.5,median(moh[,3]), paste0("S:",SRM05),
     cex=1.3, pos=3)
text(2.5,median(moh[,4]), paste0("S:",SRQ3),
     cex=1.3, pos=3)
text(3.5,median(moh[,5]), paste0("S:",SRC05),
     cex=1.3, pos=3)
text(4.5,median(moh[,6]), paste0("S:",SRSdF3),
     cex=1.3, pos=3)
text(5.5,median(moh[,7]), paste0("S:",SRSdObsI3),
     cex=1.3, pos=3)
text(6.5,median(moh[,8]), paste0("S:",SRSdObsC3),
     cex=1.3, pos=3)
text(7.5,median(moh[,9]), paste0("S:",SRCorObs),
     cex=1.3, pos=3)

axis(side = 1, at = 1:8, labels = c("M3","M/3","Q3","C/3","sdF", "sdI", "sdC", "corObs"), tck = -0.01, cex.axis = 1.5)
dev.off()

#SSB
png(filename = "figures/costalCodStrengthSSB.png",width = 700)

moh = data.frame(mohT[,2])
moh$M3 = mohM3[,2]
moh$M05 = mohM05[,2]
moh$Q3 = mohQ3[,2]
moh$C05 = mohC05[,2]
moh$SdF3 = mohSdF3[,2]
moh$SdObsI3 = mohSdObsI3[,2]
moh$SdObsC3 = mohSdObsC3[,2]
moh$CorObs = mohCorObs[,2]

library(plotrix)
p = 0.05
i = 1
par(mar=c(4,6,4,4))
plotCI(i,median(moh[,2]), li = quantile(moh[,2],p/2), ui =quantile(moh[,2],1-p/2), xlim = c(0.3,8), ylim = c(-0.7,1),
       ylab = "Mohn's rho", xlab = "Scenario", main = "Test power, SSB NCC cod", xaxt= "n", cex.axis = 1.5,cex.lab = 2.2,
       cex.main = 3, lwd = 2)
i = 2
plotCI(i,median(moh[,3]), li = quantile(moh[,3],p/2), ui =quantile(moh[,3],1-p/2), add = TRUE, lwd = 2)
i = 3
plotCI(i,median(moh[,4]), li = quantile(moh[,4],p/2), ui =quantile(moh[,4],1-p/2), add = TRUE, lwd = 2)
i = 4
plotCI(i,median(moh[,5]), li = quantile(moh[,5],p/2), ui =quantile(moh[,5],1-p/2), add = TRUE, lwd = 2)
i = 5
plotCI(i,median(moh[,6]), li = quantile(moh[,6],p/2), ui =quantile(moh[,6],1-p/2), add = TRUE, lwd = 2)
i = 6
plotCI(i,median(moh[,7]), li = quantile(moh[,7],p/2), ui =quantile(moh[,7],1-p/2), add = TRUE, lwd = 2)
i = 7
plotCI(i,median(moh[,8]), li = quantile(moh[,8],p/2), ui =quantile(moh[,8],1-p/2), add = TRUE, lwd = 2)
i = 8
plotCI(i,median(moh[,9]), li = quantile(moh[,9],p/2), ui =quantile(moh[,9],1-p/2), add = TRUE, lwd = 2)


abline(h = 0)
abline(h = quantile(moh[,1],p/2), col = 'green' ,lty = 2)
abline(h = quantile(moh[,1],1-p/2), col = 'green' ,lty = 2)

S = round(mean(moh[,1]<quantile(moh[,1],p/2) | moh[,1]>quantile(moh[,1],1-p/2)),2)
SRM3 = round(mean(moh[,2]<quantile(moh[,1],p/2) | moh[,2]>quantile(moh[,1],1-p/2)),2)
SRM05 = round(mean(moh[,3]<quantile(moh[,1],p/2) | moh[,3]>quantile(moh[,1],1-p/2)),2)
SRQ3 = round(mean(moh[,4]<quantile(moh[,1],p/2) | moh[,4]>quantile(moh[,1],1-p/2)),2)
SRC05 = round(mean(moh[,5]<quantile(moh[,1],p/2) | moh[,5]>quantile(moh[,1],1-p/2)),2)
SRSdF3 = round(mean(moh[,6]<quantile(moh[,1],p/2) | moh[,6]>quantile(moh[,1],1-p/2)),2)
SRSdObsI3 = round(mean(moh[,7]<quantile(moh[,1],p/2) | moh[,7]>quantile(moh[,1],1-p/2)),2)
SRSdObsC3 = round(mean(moh[,8]<quantile(moh[,1],p/2) | moh[,8]>quantile(moh[,1],1-p/2)),2)
SRCorObs = round(mean(moh[,9]<quantile(moh[,1],p/2) | moh[,9]>quantile(moh[,1],1-p/2)),2)

text(0.5,min(0.5,median(moh[,2])), paste0("S:",SRM3),
     cex=1.3, pos=3)
text(1.5,median(moh[,3]), paste0("S:",SRM05),
     cex=1.3, pos=3)
text(2.5,median(moh[,4]), paste0("S:",SRQ3),
     cex=1.3, pos=3)
text(3.5,median(moh[,5]), paste0("S:",SRC05),
     cex=1.3, pos=3)
text(4.5,median(moh[,6]), paste0("S:",SRSdF3),
     cex=1.3, pos=3)
text(5.5,median(moh[,7]), paste0("S:",SRSdObsI3),
     cex=1.3, pos=3)
text(6.5,median(moh[,8]), paste0("S:",SRSdObsC3),
     cex=1.3, pos=3)
text(7.5,median(moh[,9]), paste0("S:",SRCorObs),
     cex=1.3, pos=3)

axis(side = 1, at = 1:8, labels = c("M3","M/3","Q3","C/3","sdF", "sdI", "sdC", "corObs"), tck = -0.01, cex.axis = 1.5)
dev.off()

#F-bar
png(filename = "figures/costalCodStrengthF.png",width = 700)

moh = data.frame(mohT[,3])
moh$M3 = mohM3[,3]
moh$M05 = mohM05[,3]
moh$Q3 = mohQ3[,3]
moh$C05 = mohC05[,3]
moh$SdF3 = mohSdF3[,3]
moh$SdObsI3 = mohSdObsI3[,3]
moh$SdObsC3 = mohSdObsC3[,3]
moh$CorObs = mohCorObs[,3]

library(plotrix)
p = 0.05
par(mar=c(4,6,4,4))
i = 1
plotCI(i,median(moh[,2]), li = quantile(moh[,2],p/2), ui =quantile(moh[,2],1-p/2), xlim = c(0.3,8), ylim = c(-0.7,1),
       ylab = "Mohn's rho", xlab = "Scenario", main = "Test power, F-bar NCC cod", xaxt= "n", cex.axis = 1.5,cex.lab = 2.2,
       cex.main = 3, lwd = 2)
i = 2
plotCI(i,median(moh[,3]), li = quantile(moh[,3],p/2), ui =quantile(moh[,3],1-p/2), add = TRUE, lwd = 2)
i = 3
plotCI(i,median(moh[,4]), li = quantile(moh[,4],p/2), ui =quantile(moh[,4],1-p/2), add = TRUE, lwd = 2)
i = 4
plotCI(i,median(moh[,5]), li = quantile(moh[,5],p/2), ui =quantile(moh[,5],1-p/2), add = TRUE, lwd = 2)
i = 5
plotCI(i,median(moh[,6]), li = quantile(moh[,6],p/2), ui =quantile(moh[,6],1-p/2), add = TRUE, lwd = 2)
i = 6
plotCI(i,median(moh[,7]), li = quantile(moh[,7],p/2), ui =quantile(moh[,7],1-p/2), add = TRUE, lwd = 2)
i = 7
plotCI(i,median(moh[,8]), li = quantile(moh[,8],p/2), ui =quantile(moh[,8],1-p/2), add = TRUE, lwd = 2)
i = 8
plotCI(i,median(moh[,9]), li = quantile(moh[,9],p/2), ui =quantile(moh[,9],1-p/2), add = TRUE, lwd = 2)


abline(h = 0)
abline(h = quantile(moh[,1],p/2), col = 'green' ,lty = 2)
abline(h = quantile(moh[,1],1-p/2), col = 'green' ,lty = 2)

S = round(mean(moh[,1]<quantile(moh[,1],p/2) | moh[,1]>quantile(moh[,1],1-p/2)),2)
SRM3 = round(mean(moh[,2]<quantile(moh[,1],p/2) | moh[,2]>quantile(moh[,1],1-p/2)),2)
SRM05 = round(mean(moh[,3]<quantile(moh[,1],p/2) | moh[,3]>quantile(moh[,1],1-p/2)),2)
SRQ3 = round(mean(moh[,4]<quantile(moh[,1],p/2) | moh[,4]>quantile(moh[,1],1-p/2)),2)
SRC05 = round(mean(moh[,5]<quantile(moh[,1],p/2) | moh[,5]>quantile(moh[,1],1-p/2)),2)
SRSdF3 = round(mean(moh[,6]<quantile(moh[,1],p/2) | moh[,6]>quantile(moh[,1],1-p/2)),2)
SRSdObsI3 = round(mean(moh[,7]<quantile(moh[,1],p/2) | moh[,7]>quantile(moh[,1],1-p/2)),2)
SRSdObsC3 = round(mean(moh[,8]<quantile(moh[,1],p/2) | moh[,8]>quantile(moh[,1],1-p/2)),2)
SRCorObs = round(mean(moh[,9]<quantile(moh[,1],p/2) | moh[,9]>quantile(moh[,1],1-p/2)),2)


text(0.5,min(0.5,median(moh[,2])), paste0("S:",SRM3),
     cex=1.3, pos=3)
text(1.5,median(moh[,3]), paste0("S:",SRM05),
     cex=1.3, pos=3)
text(2.5,median(moh[,4]), paste0("S:",SRQ3),
     cex=1.3, pos=3)
text(3.5,median(moh[,5]), paste0("S:",SRC05),
     cex=1.3, pos=3)
text(4.5,median(moh[,6]), paste0("S:",SRSdF3),
     cex=1.3, pos=3)
text(5.5,median(moh[,7]), paste0("S:",SRSdObsI3),
     cex=1.3, pos=3)
text(6.5,median(moh[,8]), paste0("S:",SRSdObsC3),
     cex=1.3, pos=3)
text(7.5,median(moh[,9]), paste0("S:",SRCorObs),
     cex=1.3, pos=3)

axis(side = 1, at = 1:8, labels = c("M3","M/3","Q3","C/3","sdF", "sdI", "sdC", "corObs"), tck = -0.01, cex.axis = 1.5)
dev.off()







