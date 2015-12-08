#======================================================
########## Elevation Synchrony calculation ###################
#======================================================

# originally part of SPCD_MixedEffects_Sparkplots.R, but outsourced here so I can run it on all data simultaneously


require(RColorBrewer)
require(reshape)
require(dplyr)
require(MASS)
require(lme4)
require(lmerTest)
require(dplR)
require(nlme)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############# Required Functions #######################################

# core_average() & core_clean () (same as _average, but preserves all cores rather than average by tree)
source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/fn_core_average.R")
# core_trunc() = truncates cores and detrends them
source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/fn_core_trunc.R")
# error_bars()
source("/Users/leeanderegg/Desktop/R functions, general/error_bars.R")
# optimized pairs function
source("/Users/leeanderegg/Desktop/ZuurMixedModelling/AllRCode/HighstatLibV6.R")

#------------------ Function to identify chrons with internal NAs ------
internalNA <- function(x){
  x.na <- is.na(x)
  x.ok <- which(!x.na)
  n.ok <- length(x.ok)
  first.ok <- x.ok[1]
  last.ok <- x.ok[n.ok]
  if (last.ok - first.ok + 1 > n.ok) return(T)
  else { return (F)}
}

#--0 Function to identify first NA to appear in a chronology
firstNA <- function(x){
  x.na <- is.na(x)
  ifelse(length(which(x.na==T))==0, return(length(x)), return(min(which(x.na==T))))
}


##_________Function to calculate pairwise correlations__________________
paircors <- function(tmpdata, use.method="pairwise.complete.obs", summaries.only=F){
  
  ## Mid Elevation
  tmpdata.M <- tmpdata[,grep("-M-",colnames(tmpdata))]
  r.M<- cor(tmpdata.M, use = use.method) # get corr matrix
  y.M<- which(lower.tri(r.M), TRUE) # get values to subset to lower triangular
  z.M <- data.frame(Tree1 = rownames(r.M)[y.M[, 1]],Tree2 = colnames(r.M)[y.M[, 2]],cor = r.M[y.M], Elev = rep("M", times=length(r.M[y.M]))) # make a dataframe of each pairwise correlation
  # make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
  summary.M <-c(mean(z.M[,"cor"], na.rm=TRUE),sd(z.M[,"cor"], na.rm=TRUE), ncol(tmpdata.M))
  
  ## High Elevation
  tmpdata.H <- tmpdata[,grep("-H-",colnames(tmpdata))]
  r.H<- cor(tmpdata.H, use = use.method) # get corr matrix
  y.H<- which(lower.tri(r.H), TRUE) # get values to subset to lower triangular
  z.H <- data.frame(Tree1 = rownames(r.H)[y.H[, 1]],Tree2 = colnames(r.H)[y.H[, 2]],cor = r.H[y.H], Elev = rep("H", times=length(r.H[y.H]))) # make a dataframe of each pairwise correlation
  # make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
  summary.H <-c(mean(z.H[,"cor"], na.rm=TRUE),sd(z.H[,"cor"], na.rm=TRUE), ncol(tmpdata.H))
  ## Low Elevation
  if(length(grep("-L", colnames(tmpdata))>0)){
    tmpdata.L <- tmpdata[,grep("-L-",colnames(tmpdata))]
    r.L<- cor(tmpdata.L, use = use.method) # get corr matrix
    y.L<- which(lower.tri(r.L), TRUE) # get values to subset to lower triangular
    z.L <- data.frame(Tree1 = rownames(r.L)[y.L[, 1]],Tree2 = colnames(r.L)[y.L[, 2]],cor = r.L[y.L], Elev = rep("L", times=length(r.L[y.L]))) # make a dataframe of each pairwise correlation
    # make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
    summary.L <-c(mean(z.L[,"cor"], na.rm=TRUE),sd(z.L[,"cor"], na.rm=TRUE), ncol(tmpdata.L))
    sumstats <- cbind(c("L","M","H"),data.frame(rbind(summary.L,summary.M,summary.H)))
    colnames(sumstats) <- c("Elevation", "Meancorr", "SDcorr","Ntrees")
    indcors <- rbind(z.L,z.M,z.H)
  }
  else{
    sumstats <- cbind(c("M","H"),data.frame(rbind(summary.M,summary.H)))
    colnames(sumstats) <- c("Elevation", "Meancorr", "SDcorr","Ntrees")
    indcors <- rbind(z.M,z.H)
  }
  
  if(summaries.only==T){
    return(sumstats)
  }
  else{
    return(indcors)
  }
}


#_________________________________________________


#  Colorado
########## loading in the averaged ring widths (made via core_average() in Competition Analysis)
# note: I'm taking off the first column of all data because it's titled 'X'
CO_PIPOrw <- read.csv("PIPO_avgRingWidth_2-19-15.csv", header=T)[,-1]
# which(apply(CO_PIPOrw[,grep("X", colnames(CO_PIPOrw))],MARGIN = 1, FUN=internalNA)==T)
CO_POTRrw <- read.csv("POTR_avgRingWidth_2-19-15.csv", header=T)[,-1]
# which(apply(CO_POTRrw[,grep("X", colnames(CO_POTRrw))],MARGIN = 1, FUN=internalNA)==T)
CO_ABLArw <- read.csv("ABLA_avgRingWidth_2-19-15.csv", header=T)[,-1]
# which(apply(CO_ABLArw[,grep("X", colnames(CO_ABLArw))],MARGIN = 1, FUN=internalNA)==T)
#CO_ABLArw[which(CO_ABLArw$Tag=="ABLA-M-04L"),] # Pre 1950 bad, but not deleted in this dataset. deleted in one core but not the other?


# ### Traits POTR and PIPO cores from 2014 ### 
# PIPO14rw <- read.csv("PIPO_Traits_RWdata_all_v1.csv", header=T)
# POTR14rw <- read.csv("POTR_Traits_RWdata_all_v1.csv", header=T)
# #####


#+++++++ Montana ++++++++
MT_TSHErw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14//MT-TSHE-avgRingWidth_7-31-15.csv", header=T)[,-1]
# which(apply(MT_TSHErw[,grep("X", colnames(MT_TSHErw))],MARGIN = 1, FUN=internalNA)==T)
MT_PSMErw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-PSME-avgRingWidth_7-31-15.csv", header=T)[,-1]
# which(apply(MT_PSMErw[,grep("X", colnames(MT_PSMErw))],MARGIN = 1, FUN=internalNA)==T)
MT_ABLArw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-ABLA-avgRingWidth_7-31-15.csv", header=T)[,-1]
# which(apply(MT_ABLArw[,grep("X", colnames(MT_ABLArw))],MARGIN = 1, FUN=internalNA)==T)



#+++++++ Washington ++++++++
WA_TSHErw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-TSHE-avgRingWidth_7-31-15.csv", header=T)[,-1]
# gotta infill some NAs at the beginning of a core in order to caclulate BAI from the out in
WA_TSHErw[1,c("X2013","X2012","X2011","X2010","X2009")] <- mean(as.numeric(WA_TSHErw[1,paste0("X",seq(2008,1999,by=-1))]))
# which(apply(WA_TSHErw[,grep("X", colnames(WA_TSHErw))],MARGIN = 1, FUN=internalNA)==T)
WA_PSMErw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-PSME-avgRingWidth_7-31-15.csv", header=T)[,-1]
# which(apply(WA_PSMErw[,grep("X", colnames(WA_PSMErw))],MARGIN = 1, FUN=internalNA)==T)
WA_ABLArw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-ABLA-avgRingWidth_7-31-15.csv", header=T)[,-1]
# which(apply(WA_ABLArw[,grep("X", colnames(WA_ABLArw))],MARGIN = 1, FUN=internalNA)==T)
## Sweet! No internal NA sections from shitty cores. HOWEVER: prob need to check CoreTracker





### CREATING THE MOST CONSERVATIVE DATASET FOR ASSESSING synchrony
    # aka, using the length of the shortest core across all elevations for a species
CO_PIPOcutraw <- CO_PIPOrw[,1:min(apply(X = CO_PIPOrw,MARGIN = 1,FUN = firstNA))] # 74 yrs
CO_POTRcutraw <- CO_POTRrw[,1:min(apply(X = CO_POTRrw,MARGIN = 1,FUN = firstNA))] # 36 yrs would lose 2 cores < 45, 3 cores >50 (2 trees <40)
CO_ABLAcutraw <- CO_ABLArw[,1:min(apply(X = CO_ABLArw,MARGIN = 1,FUN = firstNA))] # 36 yrs would lose 3 cores < 45, 6 cores >50 (2 trees <40)

MT_TSHEcutraw <- MT_TSHErw[,1:min(apply(X = MT_TSHErw,MARGIN = 1,FUN = firstNA))] # 49yrs, only 1 tree <50 
MT_PSMEcutraw <- MT_PSMErw[,1:min(apply(X = MT_PSMErw,MARGIN = 1,FUN = firstNA))] # 51 yrs
MT_ABLAcutraw <- MT_ABLArw[,1:min(apply(X = MT_ABLArw,MARGIN = 1,FUN = firstNA))] # 45 yrs, only 2 <50


WA_TSHEcutraw <- WA_TSHErw[,1:min(apply(X = WA_TSHErw,MARGIN = 1,FUN = firstNA))] # 54 yrs
WA_PSMEcutraw <- WA_PSMErw[,1:min(apply(X = WA_PSMErw,MARGIN = 1,FUN = firstNA))] # 61 yrs
WA_ABLAcutraw <- WA_ABLArw[,1:min(apply(X = WA_ABLArw,MARGIN = 1,FUN = firstNA))] # 36 yrs, 4 trees < 45, 6 trees < 50, (2 trees < 40)

# colnames(WA_ABLArw[min(apply(X = WA_ABLArw,MARGIN = 1,FUN = firstNA))])
# 2013-(min(apply(X = WA_ABLArw,MARGIN = 1,FUN = firstNA))-6)



########### Creating Detrended RWIs ################## 
# (this duplicates code in Cd-RingWidthDetrending, but oh well)
# note, there are 5 extra columns in XXrw dataframes, so -6 from the XX_XXXdet colnumbers to get the end
# and -7 from the min col for the XX_XXXXcut colnumbers because the firstNA tells me the NA column and I want the one before


#### Parameters for core_trunc
end <- "2012"     # last date of time series
method <- "Spline" # detrending method desired. c("Spline", "ModNegExp", "Mean")
dropyrs <- 3
####


dat <- CO_PIPOrw
CO_PIPOdet <- core_trunc(workingdata=dat, first=2013-(ncol(dat)-6), last=end, detmethod=method, detrend=T, drop=dropyrs)
CO_PIPOcut <- core_trunc(workingdata=dat, first=2013+7-min(apply(X = dat,MARGIN = 1,FUN = firstNA)), last=end, detmethod=method, detrend=T, drop=0)
dat <- CO_POTRrw
CO_POTRdet <- core_trunc(workingdata=dat, first=2013-(ncol(dat)-6), last=end, detmethod=method, detrend=T, drop=dropyrs)
CO_POTRcut <- core_trunc(workingdata=dat, first=2013+7-min(apply(X = dat,MARGIN = 1,FUN = firstNA)), last=end, detmethod=method, detrend=T, drop=0)
dat <- CO_ABLArw
CO_ABLAdet <- core_trunc(workingdata=dat, first=2013-(ncol(dat)-6), last=end, detmethod=method, detrend=T, drop=dropyrs)
CO_ABLAcut <- core_trunc(workingdata=dat, first=2013+7-min(apply(X = dat,MARGIN = 1,FUN = firstNA)), last=end, detmethod=method, detrend=T, drop=0)


dat <- MT_TSHErw
MT_TSHEdet <- core_trunc(workingdata=dat, first=2013-(ncol(dat)-6), last=end, detmethod=method, detrend=T, drop=dropyrs)
MT_TSHEcut <- core_trunc(workingdata=dat, first=2013+7-min(apply(X = dat,MARGIN = 1,FUN = firstNA)), last=end, detmethod=method, detrend=T, drop=0)
dat <- MT_PSMErw
MT_PSMEdet <- core_trunc(workingdata=dat, first=2013-(ncol(dat)-6), last=end, detmethod=method, detrend=T, drop=dropyrs)
MT_PSMEcut <- core_trunc(workingdata=dat, first=2013+7-min(apply(X = dat,MARGIN = 1,FUN = firstNA)), last=end, detmethod=method, detrend=T, drop=0)
dat <- MT_ABLArw
MT_ABLAdet <- core_trunc(workingdata=dat, first=2013-(ncol(dat)-6), last=end, detmethod=method, detrend=T, drop=dropyrs)
MT_ABLAcut <- core_trunc(workingdata=dat, first=2013+7-min(apply(X = dat,MARGIN = 1,FUN = firstNA)), last=end, detmethod=method, detrend=T, drop=0)



dat <- WA_TSHErw
WA_TSHEdet <- core_trunc(workingdata=dat, first=2013-(ncol(dat)-6), last=end, detmethod=method, detrend=T, drop=dropyrs)
WA_TSHEcut <- core_trunc(workingdata=dat, first=2013+7-min(apply(X = dat,MARGIN = 1,FUN = firstNA)), last=end, detmethod=method, detrend=T, drop=0)
dat <- WA_PSMErw
WA_PSMEdet <- core_trunc(workingdata=dat, first=2013-(ncol(dat)-6), last=end, detmethod=method, detrend=T, drop=dropyrs)
WA_PSMEcut <- core_trunc(workingdata=dat, first=2013+7-min(apply(X = dat,MARGIN = 1,FUN = firstNA)), last=end, detmethod=method, detrend=T, drop=0)
dat <- WA_ABLArw
WA_ABLAdet <- core_trunc(workingdata=dat, first=2013-(ncol(dat)-6), last=end, detmethod=method, detrend=T, drop=dropyrs)
WA_ABLAcut <- core_trunc(workingdata=dat, first=2013+7-min(apply(X = dat,MARGIN = 1,FUN = firstNA)), last=end, detmethod=method, detrend=T, drop=0)





###### Calculating correlations (cut and uncut) #######

CO_PIPOcors <- paircors(tmpdata=CO_PIPOdet)
CO_PIPOcorscut <- paircors(CO_PIPOcut)
CO_PIPOsumstats <- paircors(CO_PIPOcut,summaries.only = T)

CO_POTRcors <- paircors(CO_POTRdet)
CO_POTRcorscut <- paircors(CO_POTRcut)
CO_POTRsumstats <- paircors(CO_POTRcut,summaries.only = T)

CO_ABLAcors <- paircors(CO_ABLAdet)
CO_ABLAcorscut <- paircors(CO_ABLAcut)
CO_ABLAsumstats <- paircors(CO_ABLAcut,summaries.only = T)

## HUH. this is a case where the result is highly significant with cut dataset and not with uncut!!!


MT_TSHEcors <- paircors(tmpdata=MT_TSHEdet)
MT_TSHEcorscut <- paircors(MT_TSHEcut)
MT_TSHEsumstats <- paircors(MT_TSHEcut, summaries.only = T)

MT_PSMEcors <- paircors(MT_PSMEdet)
MT_PSMEcorscut <- paircors(MT_PSMEcut) # mid elevation no longer significantly higher than L
MT_PSMEsumstats <- paircors(MT_PSMEcut, summaries.only = T)

MT_ABLAcors <- paircors(MT_ABLAdet)
MT_ABLAcorscut <- paircors(MT_ABLAcut)
MT_ABLAsumstats <- paircors(MT_ABLAcut,summaries.only = T)
## HUH. SHIT! The elevation effectis REVERSE between det and cut datasets!!!!


WA_TSHEcors <- paircors(tmpdata=WA_TSHEdet)
WA_TSHEcorscut <- paircors(WA_TSHEcut)
WA_TSHEsumstats <- paircors(WA_TSHEcut, summaries.only = T)

WA_PSMEcors <- paircors(WA_PSMEdet)
WA_PSMEcorscut <- paircors(WA_PSMEcut)
WA_PSMEsumstats <- paircors(WA_PSMEcut, summaries.only = T)


WA_ABLAcors <- paircors(WA_ABLAdet)
WA_ABLAcorscut <- paircors(WA_ABLAcut)
WA_ABLAsumstats <- paircors(WA_ABLAcut,summaries.only = T)
## HUH. SHIT! The elevation effectis don't change signs, but def change significance between det and cut datasets!!!!




#### Analyzing cut and uncut synchrony data ###########

# boxplot of CO cut vs uncut:
  # doesn't change much except ABLA, which goes from ns uncut to significant cut
quartz(width=7, height=6, title = "Colorado")
par(mfrow=c(2,3), mar=c(3,3,1,1), oma=c(0,2,2,0))
boxplot(cor~Elev, CO_PIPOcors, main="PIPO")
mtext(text = "Full", side=2, line=3)
boxplot(cor~Elev, CO_POTRcors, main="POTR")
boxplot(cor~Elev, CO_ABLAcors, main="CO_ABLA")
boxplot(cor~Elev, CO_PIPOcorscut)
mtext(text = "Cut to min core length", side=2, line=3)
boxplot(cor~Elev, CO_POTRcorscut)
boxplot(cor~Elev, CO_ABLAcorscut)

# MT boxplot
# ABLA REVERSES SIGN!
quartz(width=7, height=6, title = "Montana")
par(mfrow=c(2,3), mar=c(3,3,1,1), oma=c(0,2,2,0))
boxplot(cor~Elev, MT_TSHEcors, main="MT_TSHE")
mtext(text = "Full", side=2, line=3)
boxplot(cor~Elev, MT_PSMEcors, main="MT_PSME")
boxplot(cor~Elev, MT_ABLAcors, main="MT_ABLA")
boxplot(cor~Elev, MT_TSHEcorscut)
mtext(text = "Cut to min core length", side=2, line=3)
boxplot(cor~Elev, MT_PSMEcorscut)
boxplot(cor~Elev, MT_ABLAcorscut)

# WA boxplot:
# GAH! ABLA changes significance in cut and uncut
quartz(width=7, height=6, title = "Washington")
par(mfrow=c(2,3), mar=c(3,3,1,1), oma=c(0,2,2,0))
boxplot(cor~Elev, WA_TSHEcors, main="WA_TSHE")
mtext(text = "Full", side=2, line=3)
boxplot(cor~Elev, WA_PSMEcors, main="WA_PSME")
boxplot(cor~Elev, WA_ABLAcors, main="WA_ABLA")
boxplot(cor~Elev, WA_TSHEcorscut)
mtext(text = "Cut to min core length", side=2, line=3)
boxplot(cor~Elev, WA_PSMEcorscut)
boxplot(cor~Elev, WA_ABLAcorscut)

quartz(width=4, height=4)
#pal <- c(rgb(0,0,205,150,max=255),rgb(125,38,205,180, max=255), rgb(205,0,0,150, max=255))
pal <- brewer.pal(n = 3, name="Set1")
palette(pal)
par(mar=c(5,5,3,1))
plot(MT_ABLAcors$cor~MT_ABLAcorscut$cor, pch=16, col=MT_ABLAcors$Elev, main="MT ABLA", ylab="r of full dataset", xlab="r of cut dataset")
abline(a=0,b=1)
legend("topleft",legend = c("Low", "Mid", "High"),pch=16, col=pal)



#### actual statistical analysis ####

## Colorado ##
copiposyn1 <- gls(cor~-1 + Elev, CO_PIPOcorscut, method = "ML", weights = varIdent (form = ~1|Elev))
copiposyn0 <- gls(cor~Elev, CO_PIPOcorscut, method = "ML")
copiponull <- gls(cor~1, CO_PIPOcorscut, method = "ML")
AIC(copiponull, copiposyn0, copiposyn1)
anova(copiposyn1, copiponull) # Highly significant effect of elevation p<0.0001 DECREASING
# but looks like variance changes with elev
# tack the model estimates on to sumstats dataframe
CO_PIPOsumstats$corr.mod.mean <- as.numeric(coef(copiposyn1))
CO_PIPOsumstats$corr.mod.se <- summary(copiposyn1)$tTable[,2]
CO_PIPOsumstats$corr.mod.ci2.5 <- confint(copiposyn1)[,1]
CO_PIPOsumstats$corr.mod.ci97.5 <- confint(copiposyn1)[,2]


copotrsyn1 <- gls(cor~Elev, CO_POTRcorscut, method = "ML", weights = varIdent (form = ~1|Elev))
copotrsyn0 <- gls(cor~Elev, CO_POTRcorscut, method = "ML")
copotrnull <- gls(cor~1, CO_POTRcorscut, method = "ML")
AIC(copotrnull, copotrsyn0, copotrsyn1)
anova(copotrsyn, copotrnull) # Highly significant effect of elevation p<0.0001 DEREASING
# but variance definitely changes with elev
CO_POTRsumstats$corr.mod.mean <- as.numeric(coef(copotrsyn1))
CO_POTRsumstats$corr.mod.se <- summary(copotrsyn1)$tTable[,2]
CO_POTRsumstats$corr.mod.ci2.5 <- confint(copotrsyn1)[,1]
CO_POTRsumstats$corr.mod.ci97.5 <- confint(copotrsyn1)[,2]


coablasyn1 <- gls(cor~Elev, CO_ABLAcorscut, method = "ML", weights = varIdent (form = ~1|Elev))
coablasyn0 <- gls(cor~Elev, CO_ABLAcorscut, method = "ML")
coablanull <- gls(cor~1, CO_ABLAcorscut, method = "ML")
AIC(coablanull, coablasyn0, coablasyn1)
anova(coablasyn0, coablanull) # Highly significant effect of elevation p<0.0001 DEREASING
# variance doesn't help
# add to ABLA sumstats
CO_ABLAsumstats$corr.mod.mean <- as.numeric(coef(coablasyn0))
CO_ABLAsumstats$corr.mod.se <- summary(coablasyn0)$tTable[,2]
CO_ABLAsumstats$corr.mod.ci2.5 <- confint(coablasyn0)[,1]
CO_ABLAsumstats$corr.mod.ci97.5 <- confint(coablasyn0)[,2]



## Montana ##
mttshesyn1 <- gls(cor~-1 + Elev, MT_TSHEcorscut, method = "ML", weights = varIdent (form = ~1|Elev))
mttshesyn0 <- gls(cor~Elev, MT_TSHEcorscut, method = "ML")
mttshenull <- gls(cor~1, MT_TSHEcorscut, method = "ML")
AIC(mttshenull, mttshesyn0, mttshesyn1)
anova(mttshesyn0, mttshenull) # Highly significant effect of elevation p<0.0001 DECREASING
# but looks like variance changes with elev
# tack the model estimates on to sumstats dataframe
MT_TSHEsumstats$corr.mod.mean <- as.numeric(coef(mttshesyn1))
MT_TSHEsumstats$corr.mod.se <- summary(mttshesyn1)$tTable[,2]
MT_TSHEsumstats$corr.mod.ci2.5 <- confint(mttshesyn1)[,1]
MT_TSHEsumstats$corr.mod.ci97.5 <- confint(mttshesyn1)[,2]


mtpsmesyn1 <- gls(cor~Elev, MT_PSMEcorscut, method = "ML", weights = varIdent (form = ~1|Elev))
mtpsmesyn0 <- gls(cor~Elev, MT_PSMEcorscut, method = "ML")
mtpsmenull <- gls(cor~1, MT_PSMEcorscut, method = "ML")
AIC(mtpsmenull, mtpsmesyn0, mtpsmesyn1)
anova(mtpsmesyn0, mtpsmenull) # Highly significant effect of elevation p<0.0001 INCREASING
# but variance definitely changes with elev
MT_PSMEsumstats$corr.mod.mean <- as.numeric(coef(mtpsmesyn1))
MT_PSMEsumstats$corr.mod.se <- summary(mtpsmesyn1)$tTable[,2]
MT_PSMEsumstats$corr.mod.ci2.5 <- confint(mtpsmesyn1)[,1]
MT_PSMEsumstats$corr.mod.ci97.5 <- confint(mtpsmesyn1)[,2]


mtablasyn1 <- gls(cor~Elev, MT_ABLAcorscut, method = "ML", weights = varIdent (form = ~1|Elev))
mtablasyn0 <- gls(cor~Elev, MT_ABLAcorscut, method = "ML")
mtablanull <- gls(cor~1, MT_ABLAcorscut, method = "ML")
AIC(mtablanull, mtablasyn0, mtablasyn1)
anova(mtablasyn0, mtablanull) # Highly significant effect of elevation p<0.0001 DEREASING
# variance doesn't help
# add to ABLA sumstats
MT_ABLAsumstats$corr.mod.mean <- as.numeric(coef(mtablasyn0))
MT_ABLAsumstats$corr.mod.se <- summary(mtablasyn0)$tTable[,2]
MT_ABLAsumstats$corr.mod.ci2.5 <- confint(mtablasyn0)[,1]
MT_ABLAsumstats$corr.mod.ci97.5 <- confint(mtablasyn0)[,2]


## Washington ##
watshesyn1 <- gls(cor~-1 + Elev, WA_TSHEcorscut, method = "ML", weights = varIdent (form = ~1|Elev))
watshesyn0 <- gls(cor~Elev, WA_TSHEcorscut, method = "ML")
watshenull <- gls(cor~1, WA_TSHEcorscut, method = "ML")
AIC(watshenull, watshesyn0, watshesyn1)
anova(watshesyn0, watshenull) # Highly significant effect of elevation p<0.0001 DECREASING
# no change in variance
# tack the model estimates on to sumstats dataframe
WA_TSHEsumstats$corr.mod.mean <- as.numeric(coef(watshesyn0))
WA_TSHEsumstats$corr.mod.se <- summary(watshesyn0)$tTable[,2]
WA_TSHEsumstats$corr.mod.ci2.5 <- confint(watshesyn0)[,1]
WA_TSHEsumstats$corr.mod.ci97.5 <- confint(watshesyn0)[,2]


wapsmesyn1 <- gls(cor~Elev, WA_PSMEcorscut, method = "ML", weights = varIdent (form = ~1|Elev))
wapsmesyn0 <- gls(cor~Elev, WA_PSMEcorscut, method = "ML")
wapsmenull <- gls(cor~1, WA_PSMEcorscut, method = "ML")
AIC(wapsmenull, wapsmesyn0, wapsmesyn1)
anova(wapsmenull, wapsmesyn0) # Highly significant effect of elevation p<0.0001 HUMP SHAPED
# no change in variance
WA_PSMEsumstats$corr.mod.mean <- as.numeric(coef(wapsmesyn0))
WA_PSMEsumstats$corr.mod.se <- summary(wapsmesyn0)$tTable[,2]
WA_PSMEsumstats$corr.mod.ci2.5 <- confint(wapsmesyn0)[,1]
WA_PSMEsumstats$corr.mod.ci97.5 <- confint(wapsmesyn0)[,2]


waablasyn1 <- gls(cor~Elev, WA_ABLAcorscut, method = "ML", weights = varIdent (form = ~1|Elev))
waablasyn0 <- gls(cor~Elev, WA_ABLAcorscut, method = "ML")
waablanull <- gls(cor~1, WA_ABLAcorscut, method = "ML")
AIC(waablanull, waablasyn0, waablasyn1)
anova(waablanull, waablasyn0) # Highly significant effect of elevation p<0.0001 INCREASING
# variance doesn't help
# add to ABLA sumstats
WA_ABLAsumstats$corr.mod.mean <- as.numeric(coef(waablasyn0))
WA_ABLAsumstats$corr.mod.se <- summary(waablasyn0)$tTable[,2]
WA_ABLAsumstats$corr.mod.ci2.5 <- confint(waablasyn0)[,1]
WA_ABLAsumstats$corr.mod.ci97.5 <- confint(waablasyn0)[,2]



#### plotting code that lines up with barplots in BAI-vs-ElevComp
xvals <- c(0.6666,1.66666,2.6666)
barplot(xlim=c(0,3.33), space=0.5, width = 0.66666)









####### ***** Tree Age Analysis: preliminary **** #######
# This should be outsourced to a different code at some point.
# First pass to examine why MT_ABLA might have such strange synchrony trends cut and uncut: look at tree age
#   - can do this by looking at chronology length (lots of data, conservative and noisy age estimates)
#   - or by processing the "age" data I threw down while processing the tree cores. Some of this is ok, and 
#     some I cleaned while xdating. But some done by undergrad assistants isn't very good. 

# look at the difference in chronology lengths across elevation
lengths <- apply(MT_ABLAdet, 2, firstNA)
Elevs <- MT_ABLArw$Elev
levels(Elevs) <- list(L="L", M="M", H="H")
boxplot(lengths~Elevs)

treeage <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-ABLA-TreeAge_v1.csv", header=T)
treeage$TreeUn <- paste(treeage$Elev, treeage$Tree, sep="-")
age <- c()
elev <- c()
hitpith <- c()
treename <- unique(treeage$TreeUn)
for (i in 1:length(treename)){
  tree <- treename[i]
  tmp <- treeage[which(treeage$TreeUn==tree),]
  elev[i] <- as.character(unique(tmp$Elev))
  age[i] <- mean(tmp$nyears, na.rm=T)
  hitpith[i] <- "n"
  if(tmp$Pith[1]=="y"){
    age[i] <- tmp$nyears[1]
    hitpith[i] <- "y"
  }
  if(tmp$Pith[2]=="y"){
    age[i] <- tmp$nyears[2]
    hitpith[i] <- "y"
  }
}
elev <- factor(elev)
levels(elev) <- list(L="L", M="M", H="H")
quartz(width=3, height=3)
boxplot(age~elev, col=pal, ylab="Tree Age")
boxplot(age~elev, subset=hitpith=="y", border="red3", add=T,varwidth = T)
summary(lm(age~elev, subset=hitpith=="y"))







#++++++++++++++++++++++++++++++++++++++++++++++
##### Old plotting code, to be renovated ######
#++++++++++++++++++++++++++++++++++++++++++++++




quartz(width=5, height=5)
par(mfrow=c(3,1), mar=c(2,4,3,1), oma=c(3,0,0,0))
xlim <- c(0,1)
hist(z.H[,"cor"], main = "High Elev PIPO", xlim=xlim)
hist(z.M[,"cor"], main = "Mid Elev PIPO", xlim=xlim)
hist(z.L[,"cor"], main = "Low Elev PIPO", xlim=xlim)
mtext(text="Tree Pairwise Pearson Correlation Coefficients", side=1,outer=T)


quartz(width=5, height=5)
plot(density(corr.all[,"cor"]), type="n", ylim =c(0, max(c(density(z.H$cor)$y, density(z.M$cor)$y, density(z.L$cor)$y))+.2), main = "", xlab="", yaxs = "i")
polygon(x=c(density(z.H$cor)$x,min(density(z.H$cor)$x)), y=c(density(z.H$cor)$y,0), border="black", col =rgb(0,0,205,150,max=255))
polygon(x=c(density(z.M$cor)$x,min(density(z.M$cor)$x)), y=c(density(z.M$cor)$y,0), border="black", col= rgb(125,38,205,180, max=255))
polygon(x=c(density(z.L$cor)$x,min(density(z.L$cor)$x)), y=c(density(z.L$cor)$y,0), border="black", col=rgb(205,0,0,150, max=255) )
legend("topleft",legend=c("High", "Mid", "Low"),fill=c(rgb(0,0,205,150,max=255), rgb(125,38,205,180, max=255),rgb(205,0,0,150, max=255)), border="black",title="Elevation", bty="n")

quartz(width=5, height=5)
plot(density(corr.all[,"cor"]), type="n", ylim =c(0, max(c(density(z.H$cor)$y, density(z.M$cor)$y, density(z.L$cor)$y))+.2), main = "", xlab="", yaxs = "i")
polygon(x=c(density(z.H1$cor)$x,min(density(z.H1$cor)$x)), y=c(density(z.H1$cor)$y,0), border="black", col =rgb(0,0,205,150,max=255))
polygon(x=c(density(z.M1$cor)$x,min(density(z.M1$cor)$x)), y=c(density(z.M1$cor)$y,0), border="black", col= rgb(125,38,205,180, max=255))
polygon(x=c(density(z.L1$cor)$x,min(density(z.L1$cor)$x)), y=c(density(z.L1$cor)$y,0), border="black", col=rgb(205,0,0,150, max=255) )
legend("topleft",legend=c("High", "Mid", "Low"),fill=c(rgb(0,0,205,150,max=255), rgb(125,38,205,180, max=255),rgb(205,0,0,150, max=255)), border="black",title="Elevation", bty="n")



quartz(width=5, height=5)
par(mfrow=c(3,1), mar=c(2,4,3,1), oma=c(3,0,0,0))
xlim <- c(-0.3,1)
ylim <- c(0,110)
hist(z.H[,"cor"], main = "High Elev POTR", xlim=xlim, ylim=ylim, breaks=12, col=rgb(0,0,205,250,max=255) )#, border="blue3") # put in rgb for "blue3
hist(z.M[,"cor"], main = "Mid Elev POTR", xlim=xlim, add=T, col=rgb(125,38,205,170, max=255) )#, border="purple3") # rgb for "purple3"
hist(z.L[,"cor"], main = "Low Elev POTR", xlim=xlim, add=T, col=rgb(205,0,0,200, max=255) )#, border="red3") # rgb for "red3"
hist(z.H[,"cor"], main = "High Elev POTR", xlim=xlim, add=T, ylim=ylim, breaks=12)#, col=rgb(0,0,205,250,max=255) )#, border="blue3") # put in rgb for "blue3
dens <- 
  
  plot(density(corr.all[,"cor"]), type="n", ylim =c(0, max(c(density(z.H$cor)$y, density(z.M$cor)$y, density(z.L$cor)$y))+.2), main = "", xlab="", yaxs = "i")
polygon(x=c(density(z.H$cor)$x,min(density(z.H$cor)$x)), y=c(density(z.H$cor)$y,0), border="black", col =rgb(0,0,205,150,max=255))
polygon(x=c(density(z.M$cor)$x,min(density(z.M$cor)$x)), y=c(density(z.M$cor)$y,0), border="black", col= rgb(125,38,205,180, max=255))
polygon(x=c(density(z.L$cor)$x,min(density(z.L$cor)$x)), y=c(density(z.L$cor)$y,0), border="black", col=rgb(205,0,0,150, max=255) )
legend("topleft",legend=c("High", "Mid", "Low"),fill=c(rgb(0,0,205,150,max=255), rgb(125,38,205,180, max=255),rgb(205,0,0,150, max=255)), border="black",title="Elevation", bty="n")


mtext(text="Tree Pairwise Pearson Correlation Coefficients", side=1,outer=T)


### a breif Foray into ggplot. But turns out this will take a fair amount of tweaking
ggplot(POTR.corr.all, aes(x=cor, colour=elevs)) + geom_density() 
#  + geom_vline(data=POTR.corr, aes(xintercept=Meancorr,  colour=Elevation),
#             linetype="dashed", size=1)

"red3", "blue3", "purple3"