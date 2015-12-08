#=======================================================
#-------------------------------------------------------
#               Making detrended RWI + clim for all transects
#               Leander DL Anderegg
#                   08/31/15
#-------------------------------------------------------
#=======================================================

## load core.average() function


rm(list=ls())
#setwd("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13")
# set working directory to my scratch file in Dropbox
library(dplR)
require(reshape)
require(dplyr)
require(plyr)
require(lme4)
require(lmerTest)

################################
#         load functions   
################################
# core_average()
source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/fn_core_average.R")
# core_trunc()
source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/fn_core_trunc.R")
# error_bars()
source("/Users/leeanderegg/Desktop/R functions, general/error_bars.R")
# sparkplots()
#source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/sparkplotfn.R")
### load the core.average() function from  core_averaging_detrending.R



CO_PIPOrw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/PIPO_avgRingWidth_2-19-15.csv", header=T)
CO_PIPOrw$X <- NULL
CO_POTRrw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/POTR_avgRingWidth_2-19-15.csv", header=T)
CO_POTRrw$X <- NULL
CO_ABLArw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/ABLA_avgRingWidth_2-19-15.csv", header=T)
CO_ABLArw$X <- NULL

MT_TSHErw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-TSHE-avgRingWidth_7-20-15.csv", header=T)
MT_TSHErw$X <- NULL
MT_PSMErw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-PSME-avgRingWidth_7-31-15.csv", header=T)
MT_PSMErw$X <- NULL
MT_ABLArw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-ABLA-avgRingWidth_7-31-15.csv", header=T)
MT_ABLArw$X <- NULL

WA_TSHErw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-TSHE-avgRingWidth_7-31-15.csv", header=T)
WA_TSHErw$X <- NULL
WA_PSMErw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-PSME-avgRingWidth_7-31-15.csv", header=T)
WA_PSMErw$X <- NULL
WA_ABLArw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-ABLA-avgRingWidth_7-31-15.csv", header=T)
WA_ABLArw$X <- NULL


###########################################
#           Clean tree core Data 
###########################################

#========= subsetting to desired dates============
#-------------------------------------------------

#______________________________________
### customize these to get dates
# for PRISM data from climate WNA, first full water year = 1902, last year =2012
#workingdata <-treedata1
start <- "1902"   # first date for time series
end <- "2012"     # last date of time series
method <- "Spline" # detrending method desired. c("Spline", "ModNegExp", "Mean")
dropyrs <- 3
####

#### Making both a dataframe of detrended rwi and trend-in ring widths
#rw.det <- core_trunc(workingdata=workingdata, detmethod=method, detrend=TRUE)
#rw.nodet <- core_trunc(workingdata=workingdata, detmethod=method, detrend=FALSE)


# Colorado
CO_PIPOrwi <- core_trunc(workingdata=CO_PIPOrw, first=start, last=end, detmethod=method, detrend=TRUE, drop=dropyrs)
CO_POTRrwi <- core_trunc(workingdata=CO_POTRrw, first=start, last=end,  detmethod=method, detrend=TRUE, drop=dropyrs)
CO_ABLArwi <- core_trunc(workingdata=CO_ABLArw, first=start, last=end,  detmethod=method, detrend=TRUE, drop=dropyrs)

# Montanta
MT_TSHErwi <- core_trunc(workingdata=MT_TSHErw, first="1910", last=end,  detmethod=method, detrend=TRUE, drop=dropyrs)
### too short time series, first year = 1910
MT_PSMErwi <- core_trunc(workingdata=MT_PSMErw, first=start, last=end,  detmethod=method, detrend=TRUE, drop=dropyrs)
MT_ABLArwi <- core_trunc(workingdata=MT_ABLArw, first=start, last=end,  detmethod=method, detrend=TRUE, drop=dropyrs)

WA_TSHErwi <- core_trunc(workingdata=WA_TSHErw, first=start, last=end,  detmethod=method, detrend=TRUE, drop=dropyrs)
WA_PSMErwi <- core_trunc(workingdata=WA_PSMErw, first=start, last=end,  detmethod=method, detrend=TRUE, drop=dropyrs)
WA_ABLArwi <- core_trunc(workingdata=WA_ABLArw, first="1903", last=end,  detmethod=method, detrend=TRUE, drop=dropyrs)
## too short time series, first year =1903

#test <- data.frame(CO_PIPOrw[,1:5], t(CO_PIPOrwi))



# ## Making master chronologies
# # or using the biweighted robust mean (more the industry standard)
# PIPO.masterbw <- chron(CO_PIPOrwi, prefix="PI")
# POTR.masterbw <- chron(CO_POTRrwi, prefix="PO")
# ABLA.masterbw <- chron(CO_ABLArwi, prefix="AB")
# plot(PIstd~rownames(PIPO.masterbw), PIPO.masterbw, type="l")
# lines(POstd~rownames(POTR.masterbw), POTR.masterbw, col="blue")
# lines(ABstd~rownames(ABLA.masterbw), ABLA.masterbw, col="green")
# 
# quartz(width=5, height=5)
# par(mar=c(4,4,1,1))
# plot(POTR.masterbw$POstd~PIPO.masterbw$PIstd, pch=16, col="red", ylim=c(0.2,1.5), xlim=c(0.2,1.5), xlab="PIPO RWI", ylab="Other Species RWI")
# points(ABLA.masterbw$ABstd~PIPO.masterbw$PIstd, pch=16, col="green")
# abline(a=0, b=1)


#--------------------------------------
#==============================================





#########################################
#        Read and process climate data
#######################################


#-----------  Colorado ----
CO_PIPOclim <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/climate_CO-PIPO_09-03-15.csv", header=T) 
CO_POTRclim <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/climate_CO-POTR_09-03-15.csv", header=T)
CO_ABLAclim <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/climate_CO-ABLA_09-03-15.csv", header=T)

MT_TSHEclimall <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths//climate_MT-TSHE_09-03-15.csv", header=T)
MT_TSHEclim <- MT_TSHEclimall[-which(MT_TSHEclimall$Year %in% c(1901:1909)),] #subset down to just post 1910
MT_PSMEclim <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths//climate_MT-PSME_09-03-15.csv", header=T)
MT_ABLAclim <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths//climate_MT-ABLA_09-03-15.csv", header=T)

WA_TSHEclim <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths//climate_WA-TSHE_09-03-15.csv", header=T)
WA_PSMEclim <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths//climate_WA-PSME_09-03-15.csv", header=T)
WA_ABLAclimall <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths//climate_WA-ABLA_09-03-15.csv", header=T)
WA_ABLAclim <- WA_ABLAclimall[-which(WA_ABLAclimall$Year %in% c(1901,1902)),] #subset down to just post 1910






## Ok, things that need to happen to combine RWI and climate data
#X 1) break apart elevations for both RWI and climate
#X 2) scale climate columns (probably)
#X 3) reverse order of climate columns so they're 2012-1902
# 4) tack them together
# 5) make them long form rather than wide form
# 5) bring in any ancilary data (competition, etc, and match it to that longform data)


#### Some functions to divide up elevations and make them into lists ###### 
rwelevs_divide <- function(rwi){
  Hs <- rwi[,grep("-H-", x = colnames(rwi))]
  Ms <- rwi[,grep("-M-", x = colnames(rwi))]
  if (length(grep("-L-", x=colnames(rwi)))>1){
    Ls <- rwi[,grep("-L-", x = colnames(rwi))]
    RWIs <- list(Hs, Ms, Ls)
  }
  else {
    RWIs <- list(Hs,Ms)
  }
  return(RWIs)
}

climelevs_divide <- function (climates, scaled =TRUE){

  climatesdrop <- climates[which(climates$Year!=1901),-which(colnames(climates) %in% c("X", "ID1", "Latitude","Longitude","Elevation"))]
  
  clims <- arrange(climatesdrop, desc(Year))
  Hclims <- clims[grep("-H", x = clims$Band),]
  Mclims <- clims[grep("-M", x = clims$Band),]
  if (length(grep("-L", x=clims$Band)>1)){
    Lclims <- clims[grep("-L", x = clims$Band),]
    CLIMs <- list(Hclims, Mclims, Lclims)
  }
  else {
    CLIMs <- list(Hclims,Mclims)
  }
  if(scaled==TRUE){
    for (i in 1:length(CLIMs)){
    CLIMs[[i]] <- data.frame(CLIMs[[i]][,1:3], apply(CLIMs[[i]][,-c(1:3)], MARGIN=2, FUN=scale))
    CLIMs[[i]][,which(is.na(CLIMs[[i]][1,]))] <- rep(0, times=nrow(CLIMs[[i]]))
    }
  }
  return(CLIMs)
}



#==========================================================
###############  Put together full rwi and climate data
#=====================================================
# Put together the full data frame with rw and clim
combine_rwclim <- function(rwi,climate,Treename="Tree", omitNA = FALSE, ...) {
  require(reshape)
  rw.det <- rwelevs_divide(rwi)
  clim <- climelevs_divide(climate, ...)
  
  if(nrow(rw.det[[1]])!= nrow(clim[[1]])) {
    print("Years don't match. Shit.")
    break()
  }
  rw.climH <- data.frame(rw.det[[1]],clim[[1]])
  rw.climM <- data.frame(rw.det[[2]],clim[[2]])
  id.vars <- colnames(clim[[1]])
  measure.vars <- colnames(rw.det[[1]]) # assuming you want detrended rwi
  rw.climH.df <- melt (rw.climH, id.vars = id.vars,variable_name=Treename)
  rw.climH.df$Elev <- rep("H", times=nrow(rw.climH.df))
  rw.climM.df <- melt (rw.climM, id.vars = id.vars,variable_name=Treename)
  rw.climM.df$Elev <- rep("M", times=nrow(rw.climM.df))
  
  rw.clim.all <- rbind(rw.climH.df, rw.climM.df)
  
  if(length(rw.det)>2){
    rw.climL <- data.frame(rw.det[[3]],clim[[3]])
    rw.climL.df <- melt (rw.climL, id.vars = id.vars,variable_name=Treename)
    rw.climL.df$Elev <- rep("L", times=nrow(rw.climL.df))
    rw.clim.all <- rbind(rw.clim.all, rw.climL.df)
  }
 colnames(rw.clim.all)[which(colnames(rw.clim.all)=="value")] <- "RWI"
#   # adding in competition factor
#   ### Breaking apart H and L
#   # note, have to change whether grepping on H or L based on elevation (grep on H for L elev)
#   if(elev=="H"){
#     rw.clim.N <- rw.clim.df[grep(pattern="L",x=rw.clim.df$Tree),] 
#     rw.clim.C <- rw.clim.df[-grep(pattern="L",x=rw.clim.df$Tree),]
#     rw.clim.N$Comp <- as.factor(rep("N", times = nrow(rw.clim.N)))
#     rw.clim.C$Comp <- as.factor(rep("C", times= nrow(rw.clim.C)))
#   }
#   else {
#     rw.clim.C <- rw.clim.df[grep(pattern="H",x=rw.clim.df$Tree),] 
#     # note, have to change whether grepping on H or L based on elevation (grep on H for L elev)
#     rw.clim.N <- rw.clim.df[-grep(pattern="H",x=rw.clim.df$Tree),]
#     rw.clim.N$Comp <- as.factor(rep("N", times = nrow(rw.clim.N)))
#     rw.clim.C$Comp <- as.factor(rep("C", times= nrow(rw.clim.C)))
#   }
  #       ======= data frames to use ==========    #
  #rw.clim.df # detrended rw and climate data by tree and year and ready to go
if (omitNA==TRUE){
  rw.clim.all.narm <- na.omit(rw.clim.all) # same df but na's removed 'casue {nlme} doesn't like them
  return(rw.clim.all.narm)
}
else{
  return(rw.clim.all)
}
}
#=======================================================




CO_PIPOrwclim <- combine_rwclim (rwi = CO_PIPOrwi, climate = CO_PIPOclim)
CO_PIPOrwclim_unscaled <- combine_rwclim(rwi=CO_PIPOrwi, climate=CO_PIPOclim, scaled=FALSE)
write.csv(CO_PIPOrwclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/RW_Climdata round 2/CO_PIPO_RWIClim_101015.csv")
write.csv(CO_PIPOrwclim_unscaled, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/RW_Climdata round 2/CO_PIPO_RWI_ClimUnscaled_101015.csv")
CO_POTRrwclim <- combine_rwclim (rwi = CO_POTRrwi, climate = CO_POTRclim)
write.csv(CO_POTRrwclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/RW_Climdata round 2/CO_POTR_RWIClim_101015.csv")
CO_ABLArwclim <- combine_rwclim (rwi = CO_ABLArwi, climate = CO_ABLAclim)
write.csv(CO_ABLArwclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/RW_Climdata round 2/CO_ABLA_RWIClim_101015.csv")


MT_TSHErwclim <- combine_rwclim (rwi = MT_TSHErwi, climate = MT_TSHEclim)
write.csv(MT_TSHErwclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/RW_Climdata round 2/MT_TSHE_RWIClim_101015.csv")
MT_PSMErwclim <- combine_rwclim (rwi = MT_PSMErwi, climate = MT_PSMEclim)
write.csv(MT_PSMErwclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/RW_Climdata round 2/MT_PSME_RWIClim_101015.csv")
MT_ABLArwclim <- combine_rwclim (rwi = MT_ABLArwi, climate = MT_ABLAclim)
write.csv(MT_ABLArwclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/RW_Climdata round 2/MT_ABLA_RWIClim_101015.csv")



WA_TSHErwclim <- combine_rwclim (rwi = WA_TSHErwi, climate = WA_TSHEclim)
write.csv(WA_TSHErwclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/RW_Climdata round 2/WA_TSHE_RWIClim_101015.csv")
WA_PSMErwclim <- combine_rwclim (rwi = WA_PSMErwi, climate = WA_PSMEclim)
write.csv(WA_PSMErwclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/RW_Climdata round 2/WA_PSME_RWIClim_101015.csv")
WA_ABLArwclim <- combine_rwclim (rwi = WA_ABLArwi, climate = WA_ABLAclim)
write.csv(WA_ABLArwclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/RW_Climdata round 2/WA_ABLA_RWIClim_101015.csv")




####
dataz <- MT_TSHErwclim
# testH <- lmer (RWI ~ CMD_an + DD5_an + snow_an + PPT_an + (1+CMD_an + DD5_an + snow_an + PPT_an|Tree), data = dataz[grep("-H", x = dataz$Band),] )
#bestH <- step (testL, lsmeans.calc = FALSE)
testH <- lm (RWI ~ CMD_an + DD5_an + snow_an + PPT_an, data = dataz[grep("-H", x = dataz$Band),] )
bestH <- stepAIC(testH)
summary(bestH)
testL <- lm (RWI ~ CMD_an + DD5_an + snow_an + PPT_an, data = dataz[grep("-L", x = dataz$Band),] )
bestL <- stepAIC(testL)
summary(bestL)

testM <- lm (RWI ~ CMD_an + DD5_an + snow_an + PPT_an, data = dataz[grep("-M", x = dataz$Band),] )
bestM <- stepAIC(testM)
summary(bestM)


tmp <- WA_ABLArwclim
plot(RWI~CMD_sm , tmp, subset=Elev=="L", pch=16, cex=.5)
points(RWI~CMD_sm , tmp, subset=Elev=="H", pch=16, cex=.5,col="blue3")
abline(lm(RWI~CMD_sm , tmp, subset=Elev=="L"))
abline(lm(RWI~CMD_sm , tmp, subset=Elev=="H"), col="blue3")

plot(RWI~I(DD5_wt + DD5_sp + DD5_sm + pyDD5_at) , tmp, subset=Elev=="L", pch=16, cex=.5)
points(RWI~I(DD5_wt + DD5_sp + DD5_sm + pyDD5_at) , tmp, subset=Elev=="H", pch=16, cex=.5,col="blue3")
abline(lm(RWI~I(DD5_wt + DD5_sp + DD5_sm + pyDD5_at) , tmp, subset=Elev=="L"))
abline(lm(RWI~I(DD5_wt + DD5_sp + DD5_sm + pyDD5_at) , tmp, subset=Elev=="H"), col="blue3")


Mypairs(Z = dataz[,c("CMD_sp","CMD_sm","CMD_at","CMD_wt","CMD_an", "DD5_an","snow_an","PPT_an", "PPT_wt", "PPT_sm")])
