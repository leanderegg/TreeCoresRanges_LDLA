##############################################
######  Processing Climate Data    ###########
#### initially from climate WNA      #########
# need to convert into water year and add a new set of variables for previous year climate



#------ Functions for moving Autumn variables to the next year, and for previous year variables
clean_climate <- function(data, nIDcols=6){
  for (i in grep("_at", colnames(data))){
    data[2:nrow(data),i] <- data[1:nrow(data)-1, i]
    data[1,i] <- NA
  }
  for (i in c(grep("10", colnames(data)), grep("11", colnames(data)), grep("12", colnames(data)))){
    data[2:nrow(data),i] <- data[1:nrow(data)-1, i]
    data[1,i] <- NA
  }

  clims <- data[,(nIDcols+1):ncol(data)]
  prevyear <- rbind(rep(NA, times=ncol(clims)), clims[-nrow(clims),])
  colnames(prevyear) <- paste0("py",colnames(clims))
  data.new <- cbind(data,prevyear)
  ## calculate some derived variables
  data.new$snow_an <- with(data.new, PAS_wt + PAS_at + PAS_sp + PAS_sm )
  data.new$DD5_an <- with(data.new, DD5_01 + DD5_02 + DD5_04 + DD5_05 + DD5_06 + DD5_07 + DD5_08 + DD5_09 + DD5_10 + DD5_11 + DD5_12)
  data.new$CMD_an <- with(data.new, CMD_wt + CMD_sp + CMD_sm + CMD_at)
  data.new$PPT_an <- with(data.new, PPT_wt + PPT_sp + PPT_sm + PPT_at)
  data.new$DD5_gs <- with(data.new, DD5_05 + DD5_06 + DD5_07 + DD5_08 + DD5_09 + DD5_10)
  data.new$CMD_gs <- with(data.new, CMD05 + CMD06 + CMD07 + CMD08 + CMD09 + CMD10)
  data.new$NFFD_an <- 
  return(data.new)
}



#_________________________________________________



#### Inport Climate Data #######
# downscaled from climateWNA on 7/29/15 for all mean locations
climateall <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WholeProjectCode/climateWNA data/Timeseries/forclimateWNA_MeanLocations_8_5_15_1901-2012MSYT.csv", header=T)
climateall$Band <- paste(climateall$ID1, climateall$ID2, sep="-")
anvars <- c(231:253) # kill the annual variables that aren't in water years
climateall <- climateall[,c(1:3,ncol(climateall),4:(ncol(climateall)-1))]
climate <- climateall[,-anvars]



## split up into individual datasets
COPIPOL <- climate[which(climate$Band=="CO-PIPO-L"),]
COPIPOM <- climate[which(climate$Band=="CO-PIPO-M"),]
COPIPOH <- climate[which(climate$Band=="CO-PIPO-H"),]
COPOTRL <- climate[which(climate$Band=="CO-POTR-L"),]
COPOTRM <- climate[which(climate$Band=="CO-POTR-M"),]
COPOTRH <- climate[which(climate$Band=="CO-POTR-H"),]
COABLAL <- climate[which(climate$Band=="CO-ABLA-L"),]
COABLAM <- climate[which(climate$Band=="CO-ABLA-M"),]
COABLAH <- climate[which(climate$Band=="CO-ABLA-H"),]

MTTSHEL <- climate[which(climate$Band=="MT-TSHE-L"),]
MTTSHEM <- climate[which(climate$Band=="MT-TSHE-M"),]
MTTSHEH <- climate[which(climate$Band=="MT-TSHE-H"),]
MTPSMEL <- climate[which(climate$Band=="MT-PSME-L"),]
MTPSMEM <- climate[which(climate$Band=="MT-PSME-M"),]
MTPSMEH <- climate[which(climate$Band=="MT-PSME-H"),]
MTABLAL <- climate[which(climate$Band=="MT-ABLA-L"),]
MTABLAM <- climate[which(climate$Band=="MT-ABLA-M"),]
MTABLAH <- climate[which(climate$Band=="MT-ABLA-H"),]

WATSHEM <- climate[which(climate$Band=="WA-TSHE-M"),]
WATSHEH <- climate[which(climate$Band=="WA-TSHE-H"),]
WAPSMEL <- climate[which(climate$Band=="WA-PSME-L"),]
WAPSMEM <- climate[which(climate$Band=="WA-PSME-M"),]
WAPSMEH <- climate[which(climate$Band=="WA-PSME-H"),]
WAABLAL <- climate[which(climate$Band=="WA-ABLA-L"),]
WAABLAM <- climate[which(climate$Band=="WA-ABLA-M"),]
WAABLAH <- climate[which(climate$Band=="WA-ABLA-H"),]

allclims <- list(COPIPOL, COPIPOM,COPIPOH,COPOTRL,COPOTRM,COPOTRH,COABLAL,COABLAM,COABLAH)
allclimsnew <- lapply(X = allclims, FUN = clean_climate)
#str(allclimsnew)
# dim(allclims[[1]])
# test <- allclims[[1]]

COPIPOLclim <- allclimsnew[[1]]
COPIPOMclim <- allclimsnew[[2]]
COPIPOHclim <- allclimsnew[[3]]
COPIPOclim <- rbind(COPIPOLclim, COPIPOMclim, COPIPOHclim)
write.csv(COPIPOclim, "climate_CO-PIPO_09-03-15.csv")

COPOTRLclim <- allclimsnew[[4]]
COPOTRMclim <- allclimsnew[[5]]
COPOTRHclim <- allclimsnew[[6]]
COPOTRclim <- rbind(COPOTRLclim, COPOTRMclim, COPOTRHclim)
write.csv(COPOTRclim, "climate_CO-POTR_09-03-15.csv")
COABLALclim <- allclimsnew[[7]]
COABLAMclim <- allclimsnew[[8]]
COABLAHclim <- allclimsnew[[9]]
COABLAclim <- rbind(COABLALclim, COABLAMclim, COABLAHclim)
write.csv(COABLAclim, "climate_CO-ABLA_09-03-15.csv")





allclimsMT <- list(MTTSHEL, MTTSHEM,MTTSHEH,MTPSMEL,MTPSMEM,MTPSMEH,MTABLAL,MTABLAM,MTABLAH)
allclimsnewMT <- lapply(X = allclimsMT, FUN = clean_climate)
# str(allclimsnew)
# dim(allclims[[1]])
# test <- allclims[[1]]

MTTSHELclim <- allclimsnewMT[[1]]
MTTSHEMclim <- allclimsnewMT[[2]]
MTTSHEHclim <- allclimsnewMT[[3]]
MTTSHEclim <- rbind(MTTSHELclim, MTTSHEMclim, MTTSHEHclim)
write.csv(MTTSHEclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/climate_MT-TSHE_09-03-15.csv")
MTPSMELclim <- allclimsnewMT[[4]]
MTPSMEMclim <- allclimsnewMT[[5]]
MTPSMEHclim <- allclimsnewMT[[6]]
MTPSMEclim <- rbind(MTPSMELclim, MTPSMEMclim, MTPSMEHclim)
write.csv(MTPSMEclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/climate_MT-PSME_09-03-15.csv")
MTABLALclim <- allclimsnewMT[[7]]
MTABLAMclim <- allclimsnewMT[[8]]
MTABLAHclim <- allclimsnewMT[[9]]
MTABLAclim <- rbind(MTABLALclim, MTABLAMclim, MTABLAHclim)
write.csv(MTABLAclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/climate_MT-ABLA_09-03-15.csv")




allclimsWA <- list( WATSHEM,WATSHEH,WAPSMEL,WAPSMEM,WAPSMEH,WAABLAL,WAABLAM,WAABLAH)
allclimsnewWA <- lapply(X = allclimsWA, FUN = clean_climate)
# str(allclimsnew)
# dim(allclims[[1]])
# test <- allclims[[1]]


WATSHEMclim <- allclimsnewWA[[1]]
WATSHEHclim <- allclimsnewWA[[2]]
WATSHEclim <- rbind( WATSHEMclim, WATSHEHclim)
write.csv(WATSHEclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/climate_WA-TSHE_09-03-15.csv")
WAPSMELclim <- allclimsnewWA[[3]]
WAPSMEMclim <- allclimsnewWA[[4]]
WAPSMEHclim <- allclimsnewWA[[5]]
WAPSMEclim <- rbind(WAPSMELclim, WAPSMEMclim, WAPSMEHclim)
write.csv(WAPSMEclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/climate_WA-PSME_09-03-15.csv")
WAABLALclim <- allclimsnewWA[[6]]
WAABLAMclim <- allclimsnewWA[[7]]
WAABLAHclim <- allclimsnewWA[[8]]
WAABLAclim <- rbind(WAABLALclim, WAABLAMclim, WAABLAHclim)
write.csv(WAABLAclim, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/climate_WA-ABLA_09-03-15.csv")


