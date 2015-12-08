####### Futzing with climate data ##########
## this was done looong ago, before I redownloaded climate WNA data (8/15)
## now the climateWNA data looks better, so have never revisited this code
# But the take home seemed to be:
# 1) Lizard head isn't very long but is pretty high quality for most fields
# 2) Cortez is in F rather than C and is pretty damn messy


lizard <- read.csv("/Users/leeanderegg/Desktop/PRISM check/LizardHeadSNOTEL.csv", header=TRUE)
PIPO.master <- read.csv("MASTER_CO-PIPO_02072014.csv", header=T)


yrs <- unique(lizard$year)
precip <- c()
tmax <- c()
tmin <- c()
tavg <- c()
for(i in 1:length(yrs)){
  yr <- yrs[i]
  precip[i]<- max(lizard$prec[lizard$year==yr], na.rm=T)
  tmax[i] <- max(lizard$tmax[lizard$year==yr], na.rm=T)
  tmin[i] <- max(lizard$tmin[lizard$year==yr], na.rm=T)
  tavg[i] <- max(lizard$tavg[lizard$year==yr], na.rm=T)
}

precip[which(precip=="-Inf")] <- NA
tmax[which(tmax=="-Inf")] <- NA
tmin[which(tmin=="-Inf")] <- NA
tavg[which(tavg=="-Inf")] <- NA
prism <- climH[which(climH$Year %in% yrs),]
tr <- rev(L.master[which(names(L.master) %in% yrs)])

quartz(width=7, height=4)
par(mar=c(5,4,1,1))
plot(y=precip, x=yrs, type="l", lwd=3)
points(y=prism$PPTAn*0.039, x=prism$Year, type="l", lwd=3, col="darkred")
points(y=((tr-1)*10)+29, x=yrs[1:33], type="l", lwd=3, col="green3")
legend("topright",legend=c("PRISM", "SNOTEL", "Tree RWI"), lwd=3, col=c("darkred", "black", "green3" ))
points(y=cz$PPTAn[which(cz$years %in% yrs)], x=yrs, type="l", lwd=3, col="blue")
#require(stringr)


## Plotting master chronologies
quartz(width=7, height=4)
par(mar=c(5,4,1,1))
plot(PIPO.L~X, data=PIPO.master, type="l", lwd=2, col="red3", ylab="Detrended Master Chronologies", xlab="Year", ylim=c(0, 2.5), yaxt="n")
points(PIPO.M+0.5~X, data=PIPO.master, type="l", lwd=2, col="purple3")
points(PIPO.H+1~X, data=PIPO.master, type="l", lwd=2, col="blue3")
legend("bottom", legend=c("Low", "Mid", "High"), lwd=2, col= c("red3", "purple3", "blue3"), ncol=3, bty="n")
abline(h=1, lty=2, col="gray")
abline(h=1.5, lty=2, col="gray")
abline(h=2, lty=2, col="gray")





#didn't work
###dates <- matrix(data=rep(0,times=3), nrow=1)
#for (i in 1:nrow(lizard)){
#dates[i,1:3] <- str_split (lizard[i,1], pattern="/")
#}


## bringing in Cortez climate data
cortez <- read.csv("/Users/leeanderegg/Desktop/PRISM check/Cortez_monthly.csv", header=T, na.strings=-9999)
# get rid of rest of -999.9s
for (i in 1:ncol(cortez)){
  cortez[which(cortez[,i]==-999.9),i] <- NA
}
years <- unique(cortez$YEAR)
MaxSnow <- c()
PPTGS <- c()
PPTWin <- c()
PPTAn <- c()
TmaxAn <- c()
TminAn <- c()
TaveAn <- c()
TaveWin <- c()
TaveGS <- c()
FCconvert <- function(x){tmp <- (x-32)*5/9; return(tmp)}
for (i in 2:length(years)){ #start with 1931 to get water year
  MaxSnow[i] <- max(cortez$MXSD[which(cortez$YEAR==years[i])])
  PPTGS[i] <- 2.54/100*(sum(cortez$TPCP[which(cortez$YEAR==years[i] & cortez$MONTH %in% c(4,5,6,7,8,9))]))
  PPTWin[i] <- 2.54/100*(sum(cortez$TPCP[which(cortez$YEAR==years[i] & cortez$MONTH %in% c(1,2,3))]) + sum(cortez$TPCP[which(cortez$YEAR==years[i-1] & cortez$MONTH %in% 10:12)]))
  PPTAn[i] <- 2.54/100*(sum(cortez$TPCP[which(cortez$YEAR==years[i] & cortez$MONTH %in% 1:9)]) + sum(cortez$TPCP[which(cortez$YEAR==years[i-1] & cortez$MONTH %in% 10:12)]))
  TmaxAn[i] <- FCconvert(max(cortez$EMXT[which(cortez$YEAR==years[i])]))
  TminAn[i] <- FCconvert(max(cortez$EMNT[which(cortez$YEAR==years[i])]))
  TaveAn[i] <- FCconvert((sum(cortez$MNTM[which(cortez$YEAR==years[i] & cortez$MONTH %in% 1:9)]) + sum(cortez$MNTM[which(cortez$YEAR==years[i-1] & cortez$MONTH %in% 10:12)]))/12)
  TaveWin[i] <- FCconvert((sum(cortez$MNTM[which(cortez$YEAR==years[i] & cortez$MONTH %in% 1:3)]) + sum(cortez$MNTM[which(cortez$YEAR==years[i-1] & cortez$MONTH %in% 10:12)]))/6)
  TaveGS[i] <- FCconvert(sum(cortez$MNTM[which(cortez$YEAR==years[i] & cortez$MONTH %in% 4:9)])/6)
}
cz <- data.frame (years,PPTGS,PPTWin,PPTAn,TmaxAn,TminAn,TaveAn,TaveWin,TaveGS)  
  

