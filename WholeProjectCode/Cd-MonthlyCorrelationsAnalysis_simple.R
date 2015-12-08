
sterr <- function(x){
  se <- sd(x, na.rm=T)/sqrt(length(x[!is.na(x)]))
  return(se)
}

Band <- "L"
dataz <- CO_PIPOrwclim[which(CO_PIPOrwclim$Elev==Band),]
tempcols <- c("pyTave10", "pyTave11", "pyTave12", "pyTave01", "pyTave02", "pyTave03", "pyTave04", "pyTave05", "pyTave06", "pyTave07", "pyTave08", "pyTave09", "Tave10", "Tave11", "Tave12", "Tave01", "Tave02", "Tave03", "Tave04", "Tave05", "Tave06", "Tave07", "Tave08", "Tave09")
PPTcols <- c("pyPPT10", "pyPPT11", "pyPPT12", "pyPPT01", "pyPPT02", "pyPPT03", "pyPPT04", "pyPPT05", "pyPPT06", "pyPPT07", "pyPPT08", "pyPPT09", "PPT10", "PPT11", "PPT12", "PPT01", "PPT02", "PPT03", "PPT04", "PPT05", "PPT06", "PPT07", "PPT08", "PPT09")
CMDcols <- c("pyCMD10", "pyCMD11", "pyCMD12", "pyCMD01", "pyCMD02", "pyCMD03", "pyCMD04", "pyCMD05", "pyCMD06", "pyCMD07", "pyCMD08", "pyCMD09", "CMD10", "CMD11", "CMD12", "CMD01", "CMD02", "CMD03", "CMD04", "CMD05", "CMD06", "CMD07", "CMD08", "CMD09")
Erefcols <- c("pyEref10", "pyEref11", "pyEref12", "pyEref01", "pyEref02", "pyEref03", "pyEref04", "pyEref05", "pyEref06", "pyEref07", "pyEref08", "pyEref09", "Eref10", "Eref11", "Eref12", "Eref01", "Eref02", "Eref03", "Eref04", "Eref05", "Eref06", "Eref07", "Eref08", "Eref09")
DD5cols <- c("pyDD5_10", "pyDD5_11", "pyDD5_12", "pyDD5_01", "pyDD5_02", "pyDD5_03", "pyDD5_04", "pyDD5_05", "pyDD5_06", "pyDD5_07", "pyDD5_08", "pyDD5_09", "DD5_10", "DD5_11", "DD5_12", "DD5_01", "DD5_02", "DD5_03", "DD5_04", "DD5_05", "DD5_06", "DD5_07", "DD5_08", "DD5_09")
DD_0cols <- c("pyDD_0_10", "pyDD_0_11", "pyDD_0_12", "pyDD_0_01", "pyDD_0_02", "pyDD_0_03", "pyDD_0_04", "pyDD_0_05", "pyDD_0_06", "pyDD_0_07", "pyDD_0_08", "pyDD_0_09", "DD_0_10", "DD_0_11", "DD_0_12", "DD_0_01", "DD_0_02", "DD_0_03", "DD_0_04", "DD_0_05", "DD_0_06", "DD_0_07", "DD_0_08", "DD_0_09")
DD18cols <- c("pyDD18_10", "pyDD18_11", "pyDD18_12", "pyDD18_01", "pyDD18_02", "pyDD18_03", "pyDD18_04", "pyDD18_05", "pyDD18_06", "pyDD18_07", "pyDD18_08", "pyDD18_09", "DD18_10", "DD18_11", "DD18_12", "DD18_01", "DD18_02", "DD18_03", "DD18_04", "DD18_05", "DD18_06", "DD18_07", "DD18_08", "DD18_09")
NFFDcols <- c("pyNFFD10", "pyNFFD11", "pyNFFD12", "pyNFFD01", "pyNFFD02", "pyNFFD03", "pyNFFD04", "pyNFFD05", "pyNFFD06", "pyNFFD07", "pyNFFD08", "pyNFFD09", "NFFD10", "NFFD11", "NFFD12", "NFFD01", "NFFD02", "NFFD03", "NFFD04", "NFFD05", "NFFD06", "NFFD07", "NFFD08", "NFFD09")
RHcols <- c("pyRH10", "pyRH11", "pyRH12", "pyRH01", "pyRH02", "pyRH03", "pyRH04", "pyRH05", "pyRH06", "pyRH07", "pyRH08", "pyRH09", "RH10", "RH11", "RH12", "RH01", "RH02", "RH03", "RH04", "RH05", "RH06", "RH07", "RH08", "RH09")
PAScols <- c("pyPAS10", "pyPAS11", "pyPAS12", "pyPAS01", "pyPAS02", "pyPAS03", "pyPAS04", "pyPAS05", "pyPAS06", "pyPAS07", "pyPAS08", "pyPAS09", "PAS10", "PAS11", "PAS12", "PAS01", "PAS02", "PAS03", "PAS04", "PAS05", "PAS06", "PAS07", "PAS08", "PAS09")



corsTemp <- data.frame(matrix(rep(NA, times=24*length(unique(dataz$Tree))),ncol = 24))
colnames(corsTemp) <- tempcols
corsPPT <- data.frame(matrix(rep(NA, times=24*length(unique(dataz$Tree))),ncol = 24))
colnames(corsPPT) <- PPTcols
corsCMD <- data.frame(matrix(rep(NA, times=24*length(unique(dataz$Tree))),ncol = 24))
colnames(corsCMD) <- CMDcols
corsEref <- data.frame(matrix(rep(NA, times=24*length(unique(dataz$Tree))),ncol = 24))
corsDD5 <- data.frame(matrix(rep(NA, times=24*length(unique(dataz$Tree))),ncol = 24))
corsDD_0 <- data.frame(matrix(rep(NA, times=24*length(unique(dataz$Tree))),ncol = 24))
corsDD18 <- data.frame(matrix(rep(NA, times=24*length(unique(dataz$Tree))),ncol = 24))
corsNFFD <- data.frame(matrix(rep(NA, times=24*length(unique(dataz$Tree))),ncol = 24))
corsRH <- data.frame(matrix(rep(NA, times=24*length(unique(dataz$Tree))),ncol = 24))
corsPAS <- data.frame(matrix(rep(NA, times=24*length(unique(dataz$Tree))),ncol = 24))



for (i in 1:length(unique(dataz$Tree))){
  tree <- unique(dataz$Tree)[i]
  for(j in 1:length(tempcols)){
    corsTemp[i,j]  <- cor(x=dataz[which(dataz$Tree==tree),tempcols[j]], y=dataz[which(dataz$Tree==tree),"RWI"], use = "pairwise.complete.obs")
    corsPPT[i,j] <- cor(x=dataz[which(dataz$Tree==tree),PPTcols[j]], y=dataz[which(dataz$Tree==tree),"RWI"], use = "pairwise.complete.obs")
    corsCMD[i,j] <- cor(x=dataz[which(dataz$Tree==tree),CMDcols[j]], y=dataz[which(dataz$Tree==tree),"RWI"], use = "pairwise.complete.obs")  
    corsEref[i,j]  <- cor(x=dataz[which(dataz$Tree==tree),Erefcols[j]], y=dataz[which(dataz$Tree==tree),"RWI"], use = "pairwise.complete.obs")
    corsDD5[i,j] <- cor(x=dataz[which(dataz$Tree==tree),DD5cols[j]], y=dataz[which(dataz$Tree==tree),"RWI"], use = "pairwise.complete.obs")
    corsDD_0[i,j] <- cor(x=dataz[which(dataz$Tree==tree),DD_0cols[j]], y=dataz[which(dataz$Tree==tree),"RWI"], use = "pairwise.complete.obs")  
    corsDD18[i,j]  <- cor(x=dataz[which(dataz$Tree==tree),DD18cols[j]], y=dataz[which(dataz$Tree==tree),"RWI"], use = "pairwise.complete.obs")
    corsNFFD[i,j] <- cor(x=dataz[which(dataz$Tree==tree),NFFDcols[j]], y=dataz[which(dataz$Tree==tree),"RWI"], use = "pairwise.complete.obs")
    corsRH[i,j] <- cor(x=dataz[which(dataz$Tree==tree),RHcols[j]], y=dataz[which(dataz$Tree==tree),"RWI"], use = "pairwise.complete.obs")  
    corsPAS[i,j] <- cor(x=dataz[which(dataz$Tree==tree),PAScols[j]], y=dataz[which(dataz$Tree==tree),"RWI"], use = "pairwise.complete.obs")   
  }
}


months <- 24:1
Tempmean <- apply(corsTemp,MARGIN = 2,FUN = mean, na.rm=T)
Tempse <- apply(corsTemp, 2, sterr)
PPTmean <- apply(corsPPT,MARGIN = 2,FUN = mean, na.rm=T)
PPTse <- apply(corsPPT, 2, sterr)
CMDmean <- apply(corsCMD,MARGIN = 2,FUN = mean, na.rm=T)
CMDse <- apply(corsCMD, 2, sterr)
Erefmean <- apply(corsEref,MARGIN = 2,FUN = mean, na.rm=T)
Erefse<- apply(corsEref, 2, sterr)
DD5mean <- apply(corsDD5,MARGIN = 2,FUN = mean, na.rm=T)
DD5se<- apply(corsDD5, 2, sterr)
DD_0mean <- apply(corsDD_0,MARGIN = 2,FUN = mean, na.rm=T)
DD_0se<- apply(corsDD_0, 2, sterr)
DD18mean <- apply(corsDD18,MARGIN = 2,FUN = mean, na.rm=T)
DD18se<- apply(corsDD18, 2, sterr)
NFFDmean <- apply(corsNFFD,MARGIN = 2,FUN = mean, na.rm=T)
NFFDse<- apply(corsNFFD, 2, sterr)
RHmean <- apply(corsRH,MARGIN = 2,FUN = mean, na.rm=T)
RHse<- apply(corsRH, 2, sterr)
PASmean <- apply(corsPAS,MARGIN = 2,FUN = mean, na.rm=T)
PASse<- apply(corsPAS, 2, sterr)

# HCors <- cbind(Tempmean, Tempse, PPTmean, PPTse, CMDmean, CMDse, Erefmean, Erefse
#                 , DD5mean, DD5se, DD_0mean, DD_0se
#                 , DD18mean, DD18se, NFFDmean, NFFDse
#                 , RHmean, RHse, PASmean, PASse, months)
# MCors <- cbind(Tempmean, Tempse, PPTmean, PPTse, CMDmean, CMDse, Erefmean, Erefse
#                , DD5mean, DD5se, DD_0mean, DD_0se
#                , DD18mean, DD18se, NFFDmean, NFFDse
#                , RHmean, RHse, PASmean, PASse, months)
 LCors <- cbind(Tempmean, Tempse, PPTmean, PPTse, CMDmean, CMDse, Erefmean, Erefse
                , DD5mean, DD5se, DD_0mean, DD_0se
                , DD18mean, DD18se, NFFDmean, NFFDse
                , RHmean, RHse, PASmean, PASse, months)





## plotting all the monthly correlations
quartz(width=6, height=6)
par(mfrow=c(5,2), mar=c(2,4,0,0), oma=c(0,0,1,1), mgp=c(2,.5,0))
plot(Tempmean~months, HCors, type="l", pch=1, xaxt="n", col="red3", ylim=c(-.4,.3), ylab="Cor with Temp", xlab="Months Prior to cyNov")
abline(h=0)
error_bars("months", "Tempmean", errordata=HCors, upper="Tempse" , col="red3", length=0, lwd=2)
abline(v=9.5, lty=2)
text (x=9.5, y=.25, "CalYr", cex=.8)
abline(v=12.5, lty=3)
text( x= 12.5, y=.25, "WaterYr", cex=.8)
axis(1, at=seq(1,24,3),labels = c("Sep", "Jun", "Mar", "Dec", "pySep", "pyJun", "pyMar", "pyDec"))
lines(Tempmean~months, MCors, col="purple3")
error_bars("months", "Tempmean", errordata=MCors, upper="Tempse" , col="purple3", length = 0, lwd=2)
lines(Tempmean~months, LCors, col="green4")
error_bars("months", "Tempmean", errordata=LCors, upper="Tempse" , col="green4", length = 0, lwd=2)
legend("bottomright", legend = c("High", "Mid", "Low"), lty=1, col=c("red3", "purple3", "green4"), ncol=3, cex=.8)


plot(PPTmean~months, HCors, type="l", xaxt="n",pch=1, col="red3", ylim=c(-.2,.5), xlab="Months Prior to cyNov", ylab="Cor with PPT")
abline(h=0)
error_bars("months", "PPTmean", errordata=HCors, upper="PPTse" , col="red3", length=0, lwd=2)
abline(v=9.5, lty=2)
## text (x=9.5, y=.4, "Calendar\nYr", cex=.8)
abline(v=12.5, lty=3)
## text( x= 12.5, y=.4, "Growth/\nWater\nyr", cex=.8)
axis(1, at=seq(1,24,3),labels = c("Sep", "Jun", "Mar", "Dec", "pySep", "pyJun", "pyMar", "pyDec"))
lines(PPTmean~months, MCors, col="purple3")
error_bars("months", "PPTmean", errordata=MCors, upper="PPTse" , col="purple3", length = 0, lwd=2)
lines(PPTmean~months, LCors, col="green4")
error_bars("months", "PPTmean", errordata=LCors, upper="PPTse" , col="green4", length = 0, lwd=2)


plot(CMDmean~months, HCors, type="l", pch=1, xaxt="n", col="red3", ylim=c(-.5,.2), ylab="Cor with CMD")
abline(h=0)
error_bars("months", "CMDmean", errordata=HCors, upper="CMDse" , col="red3", length=0, lwd=2)
abline(v=9.5, lty=2)
# text (x=9.5, y=-.4, "Calendar\nYr", cex=.8)
abline(v=12.5, lty=3)
# text( x= 12.5, y=-.4, "Growth/\nWater\nyr", cex=.8)
axis(1, at=seq(1,24,3),labels = c("Sep", "Jun", "Mar", "Dec", "pySep", "pyJun", "pyMar", "pyDec"))
lines(CMDmean~months, MCors, col="purple3")
error_bars("months", "CMDmean", errordata=MCors, upper="CMDse" , col="purple3", length = 0, lwd=2)
lines(CMDmean~months, LCors, col="green4")
error_bars("months", "CMDmean", errordata=LCors, upper="CMDse" , col="green4", length = 0, lwd=2)


plot(Erefmean~months, LCors, type="l", pch=1, xaxt="n", col="green4", ylab="Cor with Eref")
abline(h=0)
error_bars("months", "Erefmean", errordata=LCors, upper="Erefse" , col="green4", length=0, lwd=2)
abline(v=9.5, lty=2)
# text (x=9.5, y=-.4, "Calendar\nYr", cex=.8)
abline(v=12.5, lty=3)
# text( x= 12.5, y=-.4, "Growth/\nWater\nyr", cex=.8)
axis(1, at=seq(1,24,3),labels = c("Sep", "Jun", "Mar", "Dec", "pySep", "pyJun", "pyMar", "pyDec"))
lines(Erefmean~months, MCors, col="purple3")
error_bars("months", "Erefmean", errordata=MCors, upper="Erefse" , col="purple3", length = 0, lwd=2)
lines(Erefmean~months, HCors, col="red3")
error_bars("months", "Erefmean", errordata=HCors, upper="Erefse" , col="red3", length = 0, lwd=2)


plot(DD5mean~months, LCors, type="l", pch=1, xaxt="n", col="green4", ylab="Cor with DD5")
abline(h=0)
error_bars("months", "DD5mean", errordata=LCors, upper="DD5se" , col="green4", length=0, lwd=2)
abline(v=9.5, lty=2)
# text (x=9.5, y=-.4, "Calendar\nYr", cex=.8)
abline(v=12.5, lty=3)
# text( x= 12.5, y=-.4, "Growth/\nWater\nyr", cex=.8)
axis(1, at=seq(1,24,3),labels = c("Sep", "Jun", "Mar", "Dec", "pySep", "pyJun", "pyMar", "pyDec"))
lines(DD5mean~months, MCors, col="purple3")
error_bars("months", "DD5mean", errordata=MCors, upper="DD5se" , col="purple3", length = 0, lwd=2)
lines(DD5mean~months, HCors, col="red3")
error_bars("months", "DD5mean", errordata=HCors, upper="DD5se" , col="red3", length = 0, lwd=2)


plot(DD_0mean~months, LCors, type="l", pch=1, xaxt="n", col="green4", ylab="Cor with DD_0")
abline(h=0)
error_bars("months", "DD_0mean", errordata=LCors, upper="DD_0se" , col="green4", length=0, lwd=2)
abline(v=9.5, lty=2)
# text (x=9.5, y=-.4, "Calendar\nYr", cex=.8)
abline(v=12.5, lty=3)
# text( x= 12.5, y=-.4, "Growth/\nWater\nyr", cex=.8)
axis(1, at=seq(1,24,3),labels = c("Sep", "Jun", "Mar", "Dec", "pySep", "pyJun", "pyMar", "pyDec"))
lines(DD_0mean~months, MCors, col="purple3")
error_bars("months", "DD_0mean", errordata=MCors, upper="DD_0se" , col="purple3", length = 0, lwd=2)
lines(DD_0mean~months, HCors, col="red3")
error_bars("months", "DD_0mean", errordata=HCors, upper="DD_0se" , col="red3", length = 0, lwd=2)




plot(DD18mean~months, LCors, type="l", pch=1, xaxt="n", col="green4", ylab="Cor with DD18")
abline(h=0)
error_bars("months", "DD18mean", errordata=LCors, upper="DD18se" , col="green4", length=0, lwd=2)
abline(v=9.5, lty=2)
# text (x=9.5, y=-.4, "Calendar\nYr", cex=.8)
abline(v=12.5, lty=3)
# text( x= 12.5, y=-.4, "Growth/\nWater\nyr", cex=.8)
axis(1, at=seq(1,24,3),labels = c("Sep", "Jun", "Mar", "Dec", "pySep", "pyJun", "pyMar", "pyDec"))
lines(DD18mean~months, MCors, col="purple3")
error_bars("months", "DD18mean", errordata=MCors, upper="DD18se" , col="purple3", length = 0, lwd=2)
lines(DD18mean~months, HCors, col="red3")
error_bars("months", "DD18mean", errordata=HCors, upper="DD18se" , col="red3", length = 0, lwd=2)



plot(NFFDmean~months, LCors, type="l", pch=1, xaxt="n", col="green4", ylab="Cor with NFFD")
abline(h=0)
error_bars("months", "NFFDmean", errordata=LCors, upper="NFFDse" , col="green4", length=0, lwd=2)
abline(v=9.5, lty=2)
# text (x=9.5, y=-.4, "Calendar\nYr", cex=.8)
abline(v=12.5, lty=3)
# text( x= 12.5, y=-.4, "Growth/\nWater\nyr", cex=.8)
axis(1, at=seq(1,24,3),labels = c("Sep", "Jun", "Mar", "Dec", "pySep", "pyJun", "pyMar", "pyDec"))
lines(NFFDmean~months, MCors, col="purple3")
error_bars("months", "NFFDmean", errordata=MCors, upper="NFFDse" , col="purple3", length = 0, lwd=2)
lines(NFFDmean~months, HCors, col="red3")
error_bars("months", "NFFDmean", errordata=HCors, upper="NFFDse" , col="red3", length = 0, lwd=2)


plot(RHmean~months, LCors, type="l", pch=1, xaxt="n", col="green4", ylab="Cor with RH")
abline(h=0)
error_bars("months", "RHmean", errordata=LCors, upper="RHse" , col="green4", length=0, lwd=2)
abline(v=9.5, lty=2)
# text (x=9.5, y=-.4, "Calendar\nYr", cex=.8)
abline(v=12.5, lty=3)
# text( x= 12.5, y=-.4, "Growth/\nWater\nyr", cex=.8)
axis(1, at=seq(1,24,3),labels = c("Sep", "Jun", "Mar", "Dec", "pySep", "pyJun", "pyMar", "pyDec"))
lines(RHmean~months, MCors, col="purple3")
error_bars("months", "RHmean", errordata=MCors, upper="RHse" , col="purple3", length = 0, lwd=2)
lines(RHmean~months, HCors, col="red3")
error_bars("months", "RHmean", errordata=HCors, upper="RHse" , col="red3", length = 0, lwd=2)


plot(PASmean~months, LCors, type="l", pch=1, xaxt="n", col="green4", ylab="Cor with PAS")
abline(h=0)
error_bars("months", "PASmean", errordata=LCors, upper="PASse" , col="green4", length=0, lwd=2)
abline(v=9.5, lty=2)
# text (x=9.5, y=-.4, "Calendar\nYr", cex=.8)
abline(v=12.5, lty=3)
# text( x= 12.5, y=-.4, "Growth/\nWater\nyr", cex=.8)
axis(1, at=seq(1,24,3),labels = c("Sep", "Jun", "Mar", "Dec", "pySep", "pyJun", "pyMar", "pyDec"))
lines(PASmean~months, MCors, col="purple3")
error_bars("months", "PASmean", errordata=MCors, upper="PASse" , col="purple3", length = 0, lwd=2)
lines(PASmean~months, HCors, col="red3")
error_bars("months", "PASmean", errordata=HCors, upper="PASse" , col="red3", length = 0, lwd=2)




### Just Precip and Temp ####
quartz(width=6, height=3)
par(mfrow=c(1,2), mar=c(3,3,0,0), oma=c(0,0,1,1), mgp=c(2,.5,0))

plot(PPTmean~months, HCors, type="l", xaxt="n",pch=1, col="red3", ylim=c(-.2,.5), xlab="Months Prior to cyOct", ylab="Cor with PPT")
abline(h=0)
error_bars("months", "PPTmean", errordata=HCors, upper="PPTse" , col="red3", length=0, lwd=2)
abline(v=9.5, lty=2)
text (x=9.5, y=.4, "Calendar\nYr", cex=.7)
abline(v=12.5, lty=3)
text( x= 12.5, y=.4, "Growth/\nWater\nyr", cex=.7)
axis(1, at=seq(1,24,3),labels = c("Sep", "Jun", "Mar", "Dec", "pySep", "pyJun", "pyMar", "pyDec"))
lines(PPTmean~months, MCors, col="purple3")
error_bars("months", "PPTmean", errordata=MCors, upper="PPTse" , col="purple3", length = 0, lwd=2)
lines(PPTmean~months, LCors, col="green4")
error_bars("months", "PPTmean", errordata=LCors, upper="PPTse" , col="green4", length = 0, lwd=2)


plot(Tempmean~months, HCors, type="l", pch=1, xaxt="n", col="red3", ylim=c(-.4,.3), ylab="Cor with Temp", xlab="Months Prior to cyOct")
abline(h=0)
error_bars("months", "Tempmean", errordata=HCors, upper="Tempse" , col="red3", length=0, lwd=2)
abline(v=9.5, lty=2)
text (x=9.5, y=.25, "Calendar\nYr", cex=.7)
abline(v=12.5, lty=3)
text( x= 12.5, y=.25, "Growth/\nWater\nyr", cex=.7)
axis(1, at=seq(1,24,3),labels = c("Sep", "Jun", "Mar", "Dec", "pySep", "pyJun", "pyMar", "pyDec"))
lines(Tempmean~months, MCors, col="purple3")
error_bars("months", "Tempmean", errordata=MCors, upper="Tempse" , col="purple3", length = 0, lwd=2)
lines(Tempmean~months, LCors, col="green4")
error_bars("months", "Tempmean", errordata=LCors, upper="Tempse" , col="green4", length = 0, lwd=2)
legend("bottomright", bty="n", legend = c("High", "Mid", "Low"), lty=1, col=c("red3", "purple3", "green4"), ncol=3, cex=.8)



