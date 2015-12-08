#########################################################
###        Basal Area Growth ~ F(elev + comp)          ###
###         For all TreeCoresRanges species            ###
#########################################################

# created 7/17/15 by LDLA
# The goal: take a bunch of the code I have now developed in
# CO_Traits/BAI_analysis + TreeCoresRanges Competition analysis
# and put it all together in one place to look at growth of
# all tree species in my 2013 dataset as f(elev + comp + DBH)



require(RColorBrewer)
require(reshape)
require(dplyr)
require(MASS)
require(lme4)
require(lmerTest)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############# Required Functions #######################################

# core_average() & core_clean () (same as _average, but preserves all cores rather than average by tree)
source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/fn_core_average.R")
# core_trunc() = truncates cores and detrends them
source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/fn_core_trunc.R")
# basal_area_inc() = converts ring widths to either BAI or radius if given DBH and bark
source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/fn_basal_area_inc.R")
# error_bars()
source("/Users/leeanderegg/Desktop/R functions, general/error_bars.R")
# optimized pairs function
source("/Users/leeanderegg/Desktop/ZuurMixedModelling/AllRCode/HighstatLibV6.R")

### Function to find year of specific radius ###
# note: this function gives you warnings if tree never got below certain size, but should just fill in with NAs
findyear <- function(radiusdata, radius){
  year <- min(which(radiusdata<=radius))
  ifelse(year==Inf, return(NA), return(year))
}

############ BAIseek() function #############
# takes a vector of years from findyear() and a matrix of BAI and calcualtes the n year mean BAI around that year for each tree.
# could also use to isolate years if vector years is all the same.
BAIseek <- function(years, BAI, nyears, lastyear){
  ## First calculate how many years to go before and after.
  after <-  round(0.5*(nyears-1) + 0.01, 0) # will put extra year on near side of timeline if nyears is even
  before <- round(0.5*(nyears-1), 0)
  actualYears <- seq(lastyear, 1800, by=-1)
  
  if (length(years)!= nrow(BAI)) print("You fucked up. years and BAI don't align")
  meanBAI <- c()
  YearReached <- c()
  for(i in 1:length(years)){
    #print(i)
    #print(years[i])
    if (is.na(years[i])){
      meanBAI[i] <- NA
      YearReached[i] <- NA
    } else{
      if(years[i]-after<1) {
        meanBAI[i] <- mean(BAI[1:nyears])
        YearReached[i] <- actualYears[years[i]] 
        print (paste("Tree too small:", names(years[i]), "row:",i))}
      else{
        meanBAI[i] <- mean(BAI[i,seq(years[i]-after, years[i]+before, by=1)])
        YearReached[i] <- actualYears[years[i]] 
      }
    }      
  } 
  tmp <- data.frame(meanBAI, YearReached)
  return(tmp)
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#========================================================================





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######################### Data #########################



#  Colorado
########## loading in the averaged ring widths (made via core_average() in Competition Analysis)
CO_PIPOrw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/PIPO_avgRingWidth_2-19-15.csv", header=T)
CO_PIPOinfo <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/PIPO_TreeInfo.csv", header=T)
levels(CO_PIPOinfo$Band) <- list(L="L", M="M", H="H")
CO_PIPOinfo$PairUn <- paste0(CO_PIPOinfo$Band,CO_PIPOinfo$Pair)
CO_PIPOinfo$ACF <- with(CO_PIPOinfo, (p_top.of.canopy. - p_bot.of.canopy.)/(p_top.of.canopy.-p_base.ref.))
CO_PIPOinfo$Height <- with(CO_PIPOinfo, (p_top.of.canopy.-p_base.ref.)*dist *0.01)
CO_PIPOinfo$BA_totLR <- with(CO_PIPOinfo, BA_totL + BA_totR)
CO_PIPOinfo <- CO_PIPOinfo[match(CO_PIPOrw$Tag,table = CO_PIPOinfo$Tag),]
# CO_PIPOinfo$Height[which(CO_PIPOinfo$Height==0)] <- NA # not necessary for PIPO
CO_PIPOinfo$DBHsc <- scale(CO_PIPOinfo$DBH)
CO_PIPOinfo$BA_totLRsc <- scale(CO_PIPOinfo$BA_totLR)
CO_PIPOinfo$N_Crsc <- scale(CO_PIPOinfo$N_Cr)
CO_PIPOinfo$in5_totsc <- scale(CO_PIPOinfo$in5_tot)
CO_PIPOinfo$ACFsc <- scale(CO_PIPOinfo$ACF)

CO_POTRrw <- read.csv("POTR_avgRingWidth_2-19-15.csv", header=T)
CO_POTRinfo <- read.csv("POTR_TreeInfo.csv", header=T)
levels(CO_POTRinfo$Band) <- list(L="L", M="M", H="H")
CO_POTRinfo$PairUn <- paste0(CO_POTRinfo$Band,CO_POTRinfo$Pair)
CO_POTRinfo$ACF <- with(CO_POTRinfo, (p_top.of.canopy. - p_bot.of.canopy.)/(p_top.of.canopy.-p_base.ref.))
CO_POTRinfo$Height <- with(CO_POTRinfo, (p_top.of.canopy.-p_base.ref.)*dist *0.01)
#CO_POTRinfo$Height[which(CO_POTRinfo$Height==0)] <- NA
CO_POTRinfo$BA_totLR <- with(CO_POTRinfo, BA_totL + BA_totR)
CO_POTRinfo <- CO_POTRinfo[match(CO_POTRrw$Tag,table = CO_POTRinfo$Tag),]
CO_POTRinfo$DBHsc <- scale(CO_POTRinfo$DBH)
CO_POTRinfo$BA_totLRsc <- scale(CO_POTRinfo$BA_totLR)
CO_POTRinfo$N_Crsc <- scale(CO_POTRinfo$N_Cr)
CO_POTRinfo$in5_totsc <- scale(CO_POTRinfo$in5_tot)
CO_POTRinfo$ACFsc <- scale(CO_POTRinfo$ACF)

CO_ABLArw <- read.csv("ABLA_avgRingWidth_2-19-15.csv", header=T)
CO_ABLAinfo <-read.csv("ABLA_TreeInfo.csv", header=T)
levels(CO_ABLAinfo$Band) <- list(L="L", M="M", H="H")
CO_ABLAinfo$PairUn <- paste0(CO_ABLAinfo$Band,CO_ABLAinfo$Pair)
CO_ABLAinfo$ACF <- with(CO_ABLAinfo, (p_top.of.canopy. - p_bot.of.canopy.)/(p_top.of.canopy.-p_base.ref.))
CO_ABLAinfo$ACF[25] <- 0.8 # no height but noted on datasheet
CO_ABLAinfo$Height <- with(CO_ABLAinfo, (p_top.of.canopy.-p_base.ref.)*dist *0.01)
# CO_ABLAinfo$Height[which(CO_ABLAinfo$Height==0)] <- NA
CO_ABLAinfo$BA_totLR <- with(CO_ABLAinfo, BA_totL + BA_totR)
CO_ABLAinfo <- CO_ABLAinfo[match(CO_ABLArw$Tag,table = CO_ABLAinfo$Tag),]
CO_ABLAinfo$Height[25] <- mean(CO_ABLAinfo$Height[which(CO_ABLAinfo$Band=="H")], na.rm=T) # fill with the mean tree height from this elevation. DBH is totally reasonable so I think this is a good choice. Puts in smack dab in middle of Height~DBH cloud
CO_ABLAinfo$DBHsc <- scale(CO_ABLAinfo$DBH)
CO_ABLAinfo$BA_totLRsc <- scale(CO_ABLAinfo$BA_totLR)
CO_ABLAinfo$N_Crsc <- scale(CO_ABLAinfo$N_Cr)
CO_ABLAinfo$in5_totsc <- scale(CO_ABLAinfo$in5_tot)
CO_ABLAinfo$ACFsc <- scale(CO_ABLAinfo$ACF)



#+++++++ Montana ++++++++
MT_allinfo <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_Checks//MT_focaltree_4-5.csv")
MT_TSHErw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14//MT-TSHE-avgRingWidth_7-20-15.csv", header=T)
MT_TSHEinfo <- MT_allinfo[which(MT_allinfo$Species=="TSHE"),]
levels(MT_TSHEinfo$Band) <- list(L="L", M="M", H="H")
MT_TSHEinfo$Tag <- factor(MT_TSHEinfo$Tag) #get rid of unused levels
MT_TSHEinfo$PairUn <- paste0(MT_TSHEinfo$Band,MT_TSHEinfo$Pair)
MT_TSHEinfo$ACF <- with(MT_TSHEinfo, (p_top.of.canopy. - p_bot.of.canopy.)/(p_top.of.canopy.-p_base.ref.))
MT_TSHEinfo$Height <- with(MT_TSHEinfo, (p_top.of.canopy.-p_base.ref.)*dist *0.01)
MT_TSHEinfo <- MT_TSHEinfo[match(MT_TSHErw$Tag, MT_TSHEinfo$Tag),]
MT_TSHEinfo$BA_totLR <- with(MT_TSHEinfo, BA_totL + BA_totR)
#MT_TSHEinfo$Tag[which(is.na(MT_TSHEinfo$N_Cr))] #[51] TSHE-L-13C
MT_TSHEinfo$N_Cr[51] <- 5 # based on the one fisheye of this tree. other fisheye missing
# also had to switch the bot.canopy and top.canopy for [39] TSHE-L-07C, (top=78, bot=34). did this in the MT_focaltree_4-5.csv
MT_TSHEinfo$DBHsc <- scale(MT_TSHEinfo$DBH)
MT_TSHEinfo$BA_totLRsc <- scale(MT_TSHEinfo$BA_totLR)
MT_TSHEinfo$N_Crsc <- scale(MT_TSHEinfo$N_Cr)
MT_TSHEinfo$in5_totsc <- scale(MT_TSHEinfo$in5_tot)
MT_TSHEinfo$ACFsc <- scale(MT_TSHEinfo$ACF)

MT_PSMErw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-PSME-avgRingWidth_7-20-15.csv", header=T)
MT_PSMEinfo <- MT_allinfo[which(MT_allinfo$Species=="PSME"),]
levels(MT_PSMEinfo$Band) <- list(L="L", M="M", H="H")
MT_PSMEinfo$Tag <- factor(MT_PSMEinfo$Tag) # get rid of unused levels
MT_PSMEinfo$PairUn <- paste0(MT_PSMEinfo$Band,MT_PSMEinfo$Pair)
MT_PSMEinfo$ACF <- with(MT_PSMEinfo, (p_top.of.canopy. - p_bot.of.canopy.)/(p_top.of.canopy.-p_base.ref.))
MT_PSMEinfo$Height <- with(MT_PSMEinfo, (p_top.of.canopy.-p_base.ref.)*dist *0.01)
MT_PSMEinfo <- MT_PSMEinfo[match(MT_PSMErw$Tag, MT_PSMEinfo$Tag),]
MT_PSMEinfo$BA_totLR <- with(MT_PSMEinfo, BA_totL + BA_totR)
MT_PSMEinfo$DBHsc <- scale(MT_PSMEinfo$DBH)
MT_PSMEinfo$BA_totLRsc <- scale(MT_PSMEinfo$BA_totLR)
MT_PSMEinfo$N_Crsc <- scale(MT_PSMEinfo$N_Cr)
MT_PSMEinfo$in5_totsc <- scale(MT_PSMEinfo$in5_tot)
MT_PSMEinfo$ACFsc <- scale(MT_PSMEinfo$ACF)

MT_ABLArw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-ABLA-avgRingWidth_7-20-15.csv", header=T)
MT_ABLAinfo <-MT_allinfo[which(MT_allinfo$Species=="ABLA"),]
levels(MT_ABLAinfo$Band) <- list(L="L", M="M", H="H")
MT_ABLAinfo$Tag <- factor(MT_ABLAinfo$Tag)
MT_ABLAinfo$PairUn <- paste0(MT_ABLAinfo$Band,MT_ABLAinfo$Pair)
MT_ABLAinfo$ACF <- with(MT_ABLAinfo, (p_top.of.canopy. - p_bot.of.canopy.)/(p_top.of.canopy.-p_base.ref.))
MT_ABLAinfo$Height <- with(MT_ABLAinfo, (p_top.of.canopy.-p_base.ref.)*dist *0.01)
# MT_ABLAinfo$Height[which(MT_ABLAinfo$Height==0)] <- NA
MT_ABLAinfo <- MT_ABLAinfo[match(MT_ABLArw$Tag, MT_ABLAinfo$Tag),]
MT_ABLAinfo$BA_totLR <- with(MT_ABLAinfo, BA_totL + BA_totR)
MT_ABLAinfo$DBHsc <- scale(MT_ABLAinfo$DBH)
MT_ABLAinfo$BA_totLRsc <- scale(MT_ABLAinfo$BA_totLR)
MT_ABLAinfo$N_Crsc <- scale(MT_ABLAinfo$N_Cr)
MT_ABLAinfo$in5_totsc <- scale(MT_ABLAinfo$in5_tot)
MT_ABLAinfo$ACFsc <- scale(MT_ABLAinfo$ACF)


#+++++++ Washington ++++++++
WA_allinfo <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_Checks/WA_focaltree_4-18.csv")
WA_TSHErw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-TSHE-avgRingWidth_7-20-15.csv", header=T)
WA_TSHEinfo <- WA_allinfo[which(WA_allinfo$Species=="TSHE"),]
levels(WA_TSHEinfo$Band) <- list(L="L", M="M", H="H")
WA_TSHEinfo$Tag <- factor(WA_TSHEinfo$Tag) #get rid of unused levels
WA_TSHEinfo$PairUn <- paste0(WA_TSHEinfo$Band,WA_TSHEinfo$Pair)
WA_TSHEinfo$ACF <- with(WA_TSHEinfo, (p_top.of.canopy. - p_bot.of.canopy.)/(p_top.of.canopy.-p_base.ref.))
WA_TSHEinfo$Height <- with(WA_TSHEinfo, (p_top.of.canopy.-p_base.ref.)*dist *0.01)
WA_TSHEinfo <- WA_TSHEinfo[match(WA_TSHErw$Tag, WA_TSHEinfo$Tag),]
# gotta infill some NAs at the beginning of a core in order to caclulate BAI from the out in
WA_TSHErw[1,c("X2013","X2012","X2011","X2010","X2009")] <- mean(as.numeric(WA_TSHErw[1,paste0("X",seq(2008,1999,by=-1))]))
WA_TSHEinfo$BA_totLR <- with(WA_TSHEinfo, BA_totL + BA_totR)
WA_TSHEinfo$DBHsc <- scale(WA_TSHEinfo$DBH)
WA_TSHEinfo$BA_totLRsc <- scale(WA_TSHEinfo$BA_totLR)
WA_TSHEinfo$N_Crsc <- scale(WA_TSHEinfo$N_Cr)
WA_TSHEinfo$in5_totsc <- scale(WA_TSHEinfo$in5_tot)
WA_TSHEinfo$ACFsc <- scale(WA_TSHEinfo$ACF)

WA_PSMErw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-PSME-avgRingWidth_7-20-15.csv", header=T)
WA_PSMEinfo <- WA_allinfo[which(WA_allinfo$Species=="PSME"),]
levels(WA_PSMEinfo$Band) <- list(L="L", M="M", H="H")
WA_PSMEinfo$Tag <- factor(WA_PSMEinfo$Tag) # get rid of unused levels
WA_PSMEinfo$PairUn <- paste0(WA_PSMEinfo$Band,WA_PSMEinfo$Pair)
WA_PSMEinfo$ACF <- with(WA_PSMEinfo, (p_top.of.canopy. - p_bot.of.canopy.)/(p_top.of.canopy.-p_base.ref.))
WA_PSMEinfo$Height <- with(WA_PSMEinfo, (p_top.of.canopy.-p_base.ref.)*dist *0.01)
WA_PSMEinfo <- WA_PSMEinfo[match(WA_PSMErw$Tag, WA_PSMEinfo$Tag),]
WA_PSMEinfo$BA_totLR <- with(WA_PSMEinfo, BA_totL + BA_totR)
WA_PSMEinfo$DBHsc <- scale(WA_PSMEinfo$DBH)
WA_PSMEinfo$BA_totLRsc <- scale(WA_PSMEinfo$BA_totLR)
WA_PSMEinfo$N_Crsc <- scale(WA_PSMEinfo$N_Cr)
WA_PSMEinfo$in5_totsc <- scale(WA_PSMEinfo$in5_tot)
WA_PSMEinfo$ACFsc <- scale(WA_PSMEinfo$ACF)

WA_ABLArw <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-ABLA-avgRingWidth_7-20-15.csv", header=T)
WA_ABLAinfo <-WA_allinfo[which(WA_allinfo$Species=="ABLA"),]
levels(WA_ABLAinfo$Band) <- list(L="L", M="M", H="H")
WA_ABLAinfo$Tag <- factor(WA_ABLAinfo$Tag)
WA_ABLAinfo$PairUn <- paste0(WA_ABLAinfo$Band,WA_ABLAinfo$Pair)
WA_ABLAinfo$ACF <- with(WA_ABLAinfo, (p_top.of.canopy. - p_bot.of.canopy.)/(p_top.of.canopy.-p_base.ref.))
WA_ABLAinfo$Height <- with(WA_ABLAinfo, (p_top.of.canopy.-p_base.ref.)*dist *0.01)
# WA_ABLAinfo$Height[which(WA_ABLAinfo$Height==0)] <- NA
WA_ABLAinfo <- WA_ABLAinfo[match(WA_ABLArw$Tag, WA_ABLAinfo$Tag),]
WA_ABLAinfo$BA_totLR <- with(WA_ABLAinfo, BA_totL + BA_totR)
WA_ABLAinfo$DBHsc <- scale(WA_ABLAinfo$DBH)
WA_ABLAinfo$BA_totLRsc <- scale(WA_ABLAinfo$BA_totLR)
WA_ABLAinfo$N_Crsc <- scale(WA_ABLAinfo$N_Cr)
WA_ABLAinfo$in5_totsc <- scale(WA_ABLAinfo$in5_tot)
WA_ABLAinfo$ACFsc <- scale(WA_ABLAinfo$ACF)


##################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 
# # eyeballing for outliers
# dataz <- CO_PIPOinfo
# for(i in 1:length(comps))
# {hist(dataz[,comps[i]], main=paste(dataz$State[1],dataz$Species[1],comps[i]), xlab=comps[i], breaks=10)}
# plot(Height~DBH, pch=16, col= Band, data=dataz)
# apply(dataz[,comps], 2, min)
# apply(dataz[,comps], 2, max)
# #WA_TSHEinfo$Tag[which(WA_TSHEinfo$in5_tot>12)]







#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############## Creating BAI ###################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# the approach: get the same 10yr section for each tree to take avg BAI, sd BAI, etc.

#+++++++++
## Caclulating some BAI stats from 2003-2012
interval <- 2:11 # these correspond to 2003-2012
#+++++++++


#------ Colorado ---
CO_pipobai <- basal_area_inc(ringwidths=CO_PIPOrw, DBHbark=CO_PIPOinfo,calBAI = T, DBHcol =which(colnames(CO_PIPOinfo)=='DBH'), barkcol=which(colnames(CO_PIPOinfo)=='Bark') )
#min(CO_pipobai, na.rm=T) # that seemed to do it. got a tree that's 8mm negative, but will just roll with it
CO_PIPOinfo$BAI10yr <- apply(X=CO_pipobai[,interval],MARGIN=1,FUN=mean) 
CO_PIPOinfo$BAI10yrtot <- apply (X=CO_pipobai[,interval],MARGIN=1,FUN=sum)
CO_PIPOinfo$BAI10yrmax <- apply (X=CO_pipobai[,interval],MARGIN=1,FUN=max)
CO_PIPOinfo$BAI10yrmin <- apply (X=CO_pipobai[,interval],MARGIN=1,FUN=min)
CO_PIPOinfo$BAI10yrsd <- apply (X=CO_pipobai[,interval],MARGIN=1,FUN=sd)
CO_PIPOinfo$BAI10yrCV <- with(CO_PIPOinfo, BAI10yrsd / BAI10yr)

CO_potrbai <- basal_area_inc(ringwidths=CO_POTRrw, DBHbark=CO_POTRinfo,calBAI = T, DBHcol =which(colnames(CO_POTRinfo)=='DBH'), barkcol=which(colnames(CO_POTRinfo)=='Bark') )
#min(CO_potrbai, na.rm=T) # got a tree that's 16.3mm negative. But I think that is not too bad
CO_POTRinfo$BAI10yr <- apply(X=CO_potrbai[,interval],MARGIN=1,FUN=mean) 
CO_POTRinfo$BAI10yrtot <- apply (X=CO_potrbai[,interval],MARGIN=1,FUN=sum)
CO_POTRinfo$BAI10yrmax <- apply (X=CO_potrbai[,interval],MARGIN=1,FUN=max)
CO_POTRinfo$BAI10yrmin <- apply (X=CO_potrbai[,interval],MARGIN=1,FUN=min)
CO_POTRinfo$BAI10yrsd <- apply (X=CO_potrbai[,interval],MARGIN=1,FUN=sd)
CO_POTRinfo$BAI10yrCV <- with(CO_POTRinfo, BAI10yrsd / BAI10yr)

CO_ablabai <- basal_area_inc(ringwidths=CO_ABLArw, DBHbark=CO_ABLAinfo,calBAI = T, DBHcol =which(colnames(CO_ABLAinfo)=='DBH'), barkcol=which(colnames(CO_ABLAinfo)=='Bark') )
CO_ablarad <- basal_area_inc(ringwidths=CO_ABLArw, DBHbark=CO_ABLAinfo,calBAI = F, DBHcol =which(colnames(CO_ABLAinfo)=='DBH'), barkcol=which(colnames(CO_ABLAinfo)=='Bark') )
# yeeesh. got a tree that's 54mm off.
#****** need to remove a couple trees from dataset, because SOMETHING is clearly off
# use findyear() to pull out those trees that have radius < -10 and get rid of them
CO_ablabai <- CO_ablabai[-which(apply(CO_ablarad, MARGIN=1,FUN=findyear,radius=-20)>1),]
CO_ABLAinfo$BAI10yr <- apply(X=CO_ablabai[,interval],MARGIN=1,FUN=mean)[match(CO_ABLAinfo$Tag,rownames(CO_ablabai))] 
CO_ABLAinfo$BAI10yrtot <- apply (X=CO_ablabai[,interval],MARGIN=1,FUN=sum)[match(CO_ABLAinfo$Tag,rownames(CO_ablabai))]
CO_ABLAinfo$BAI10yrmax <- apply (X=CO_ablabai[,interval],MARGIN=1,FUN=max)[match(CO_ABLAinfo$Tag,rownames(CO_ablabai))]
CO_ABLAinfo$BAI10yrmin <- apply (X=CO_ablabai[,interval],MARGIN=1,FUN=min)[match(CO_ABLAinfo$Tag,rownames(CO_ablabai))]
CO_ABLAinfo$BAI10yrsd <- apply (X=CO_ablabai[,interval],MARGIN=1,FUN=sd)[match(CO_ABLAinfo$Tag,rownames(CO_ablabai))]
CO_ABLAinfo$BAI10yrCV <- with(CO_ABLAinfo, BAI10yrsd / BAI10yr)





#------- Montana ---
MT_TSHEbai <- basal_area_inc(ringwidths=MT_TSHErw, DBHbark=MT_TSHEinfo,calBAI = T, DBHcol =which(colnames(MT_TSHEinfo)=='DBH'), barkcol=which(colnames(MT_TSHEinfo)=='Bark') )
#min(CO_pipobai, na.rm=T) # that seemed to do it. got a tree that's 3.8mm negative, but will just roll with it
MT_TSHEinfo$BAI10yr <- apply(X=MT_TSHEbai[,interval],MARGIN=1,FUN=mean) 
MT_TSHEinfo$BAI10yrtot <- apply (X=MT_TSHEbai[,interval],MARGIN=1,FUN=sum)
MT_TSHEinfo$BAI10yrmax <- apply (X=MT_TSHEbai[,interval],MARGIN=1,FUN=max)
MT_TSHEinfo$BAI10yrmin <- apply (X=MT_TSHEbai[,interval],MARGIN=1,FUN=min)
MT_TSHEinfo$BAI10yrsd <- apply (X=MT_TSHEbai[,interval], MARGIN=1, FUN=sd)
MT_TSHEinfo$BAI10yrCV <- with(MT_TSHEinfo, BAI10yrsd / BAI10yr)

MT_PSMEbai <- basal_area_inc(ringwidths=MT_PSMErw, DBHbark=MT_PSMEinfo,calBAI = T, DBHcol =which(colnames(MT_PSMEinfo)=='DBH'), barkcol=which(colnames(MT_PSMEinfo)=='Bark') )
#min(CO_pipobai, na.rm=T) # that seemed to do it. got a tree that's 1.2mm negative, but will just roll with it
MT_PSMEinfo$BAI10yr <- apply(X=MT_PSMEbai[,interval],MARGIN=1,FUN=mean) 
MT_PSMEinfo$BAI10yrtot <- apply (X=MT_PSMEbai[,interval],MARGIN=1,FUN=sum)
MT_PSMEinfo$BAI10yrmax <- apply (X=MT_PSMEbai[,interval],MARGIN=1,FUN=max)
MT_PSMEinfo$BAI10yrmin <- apply (X=MT_PSMEbai[,interval],MARGIN=1,FUN=min)
MT_PSMEinfo$BAI10yrsd <- apply (X=MT_PSMEbai[,interval], MARGIN=1, FUN=sd)
MT_PSMEinfo$BAI10yrCV <- with(MT_PSMEinfo, BAI10yrsd / BAI10yr)

MT_ABLAbai <- basal_area_inc(ringwidths=MT_ABLArw, DBHbark=MT_ABLAinfo,calBAI = T, DBHcol =which(colnames(MT_ABLAinfo)=='DBH'), barkcol=which(colnames(MT_ABLAinfo)=='Bark') )
#min(CO_pipobai, na.rm=T) # that seemed to do it. got a tree that's 4.6mm negative, but will just roll with it
MT_ABLAinfo$BAI10yr <- apply(X=MT_ABLAbai[,interval],MARGIN=1,FUN=mean) 
MT_ABLAinfo$BAI10yrtot <- apply (X=MT_ABLAbai[,interval],MARGIN=1,FUN=sum)
MT_ABLAinfo$BAI10yrmax <- apply (X=MT_ABLAbai[,interval],MARGIN=1,FUN=max)
MT_ABLAinfo$BAI10yrmin <- apply (X=MT_ABLAbai[,interval],MARGIN=1,FUN=min)
MT_ABLAinfo$BAI10yrsd <- apply (X=MT_ABLAbai[,interval], MARGIN=1, FUN=sd)
MT_ABLAinfo$BAI10yrCV <- with(MT_ABLAinfo, BAI10yrsd / BAI10yr)



#------- Washington ---
WA_TSHEbai <- basal_area_inc(ringwidths=WA_TSHErw, DBHbark=WA_TSHEinfo,calBAI = T, DBHcol =which(colnames(WA_TSHEinfo)=='DBH'), barkcol=which(colnames(WA_TSHEinfo)=='Bark') )
#min(CO_pipobai, na.rm=T) # that seemed to do it. got a tree that's 1.2mm negative, but will just roll with it
WA_TSHEinfo$BAI10yr <- apply(X=WA_TSHEbai[,interval],MARGIN=1,FUN=mean) 
WA_TSHEinfo$BAI10yrtot <- apply (X=WA_TSHEbai[,interval],MARGIN=1,FUN=sum)
WA_TSHEinfo$BAI10yrmax <- apply (X=WA_TSHEbai[,interval],MARGIN=1,FUN=max)
WA_TSHEinfo$BAI10yrmin <- apply (X=WA_TSHEbai[,interval],MARGIN=1,FUN=min)
WA_TSHEinfo$BAI10yrsd <- apply (X=WA_TSHEbai[,interval], MARGIN=1, FUN=sd)
WA_TSHEinfo$BAI10yrCV <- with(WA_TSHEinfo, BAI10yrsd / BAI10yr)

WA_PSMEbai <- basal_area_inc(ringwidths=WA_PSMErw, DBHbark=WA_PSMEinfo,calBAI = T, DBHcol =which(colnames(WA_PSMEinfo)=='DBH'), barkcol=which(colnames(WA_PSMEinfo)=='Bark') )
#min(CO_pipobai, na.rm=T) # that seemed to do it. got a tree that's 3.3mm negative, but will just roll with it
WA_PSMEinfo$BAI10yr <- apply(X=WA_PSMEbai[,interval],MARGIN=1,FUN=mean) 
WA_PSMEinfo$BAI10yrtot <- apply (X=WA_PSMEbai[,interval],MARGIN=1,FUN=sum)
WA_PSMEinfo$BAI10yrmax <- apply (X=WA_PSMEbai[,interval],MARGIN=1,FUN=max)
WA_PSMEinfo$BAI10yrmin <- apply (X=WA_PSMEbai[,interval],MARGIN=1,FUN=min)
WA_PSMEinfo$BAI10yrsd <- apply (X=WA_PSMEbai[,interval], MARGIN=1, FUN=sd)
WA_PSMEinfo$BAI10yrCV <- with(WA_PSMEinfo, BAI10yrsd / BAI10yr)

WA_ABLAbai <- basal_area_inc(ringwidths=WA_ABLArw, DBHbark=WA_ABLAinfo,calBAI = T, DBHcol =which(colnames(WA_ABLAinfo)=='DBH'), barkcol=which(colnames(WA_ABLAinfo)=='Bark') )
#min(CO_pipobai, na.rm=T) # damn. have a tree that's -11mm
WA_ABLArad <- basal_area_inc(ringwidths=WA_ABLArw, DBHbark=WA_ABLAinfo,calBAI = F, DBHcol =which(colnames(WA_ABLAinfo)=='DBH'), barkcol=which(colnames(WA_ABLAinfo)=='Bark') )
# which(apply(WA_ABLArad, MARGIN=1,FUN=findyear,radius=-10)>1)
# ABLA-M-6N is the offending party. But I think I'm going to leave him in for the time being
WA_ABLAinfo$BAI10yr <- apply(X=WA_ABLAbai[,interval],MARGIN=1,FUN=mean) 
WA_ABLAinfo$BAI10yrtot <- apply (X=WA_ABLAbai[,interval],MARGIN=1,FUN=sum)
WA_ABLAinfo$BAI10yrmax <- apply (X=WA_ABLAbai[,interval],MARGIN=1,FUN=max)
WA_ABLAinfo$BAI10yrmin <- apply (X=WA_ABLAbai[,interval],MARGIN=1,FUN=min)
WA_ABLAinfo$BAI10yrsd <- apply (X=WA_ABLAbai[,interval], MARGIN=1, FUN=sd)
WA_ABLAinfo$BAI10yrCV <- with(WA_ABLAinfo, BAI10yrsd / BAI10yr)




#write.csv(CO_PIPOinfo,"/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/"); write.csv(CO_POTRinfo,); write.csv(CO_ABLAinfo,)
#write.csv(MT_THSEinfo,); write.csv(MT_PSMEinfo,); write.csv(MT_ABLAinfo,)




####### Colinearity Check for Competative Indeces ##############
# I got concerned of corr >0.7 or if VIF>3. Suprisingly infrequent...

#  Compcols <- c( "N_Cr","BA_tot", "BA_totL", "BA_same_all","BA_sameL", "BA_other_all", "BA_otherL",
#                "DBH","in5_tot", "in5_L","in5_same","in5_sameL", "in5_other", "in5_otherL", "ACF")
#  comps <- c("N_Cr","BA_totLR","BA_sameLR","BA_otherLR", "DBH", "in5_tot","in5_same","in5_other","ACF", "Height")
#  
#  Mypairs(Z = CO_PIPOinfo[,comps])#, color=c("red3","purple","blue3")[unclass(CO_ABLAinfo$Band)])
#  corvif(dataz= CO_PIPOinfo[,c("N_Cr","BA_totLR","in5_tot", "ACF")]) 
#  # N_Cr and in5_tot real correlated
#  # ACF and BA_tot and in5_same real correlated
#  Mypairs(Z = CO_POTRinfo[,comps])#, color=c("red3","purple","blue3")[unclass(CO_ABLAinfo$Band)])
#  corvif(dataz= CO_POTRinfo[,c("N_Cr","BA_totLR","in5_tot", "ACF")]) 
#  # BA_tot & in5_tot 0.66 + N_Cr & BA_tot 0.54, but VIFs not so bad?
#  Mypairs(Z = CO_ABLAinfo[,comps])#, color=c("red3","purple","blue3")[unclass(CO_ABLAinfo$Band)])
#  corvif(dataz= CO_ABLAinfo[,c("N_Cr","BA_totLR","in5_tot", "ACF")]) 
#  # almsot all of them corr>0.6 (N_Cr,BA_tot & in5_tot all close to .7), but VIFs ok?
#  Mypairs(Z = MT_TSHEinfo[,comps])#, color=c("red3","purple","blue3")[unclass(CO_ABLAinfo$Band)])
#  corvif(dataz= MT_TSHEinfo[,c("N_Cr","BA_totLR","in5_tot", "ACF")]) 
#  # only N_Cr & in5_tot > 0.6, but VIFs all look good
#  corvif(dataz= MT_PSMEinfo[,c("N_Cr","BA_totLR","in5_tot", "ACF")]) 
#  Mypairs(Z = MT_PSMEinfo[,comps])#, color=c("red3","purple","blue3")[unclass(CO_ABLAinfo$Band)])
#  corvif(dataz= MT_PSMEinfo[,c("N_Cr","BA_totLR","in5_tot", "ACF")]) 
#  # N_Cr & in5_tot>0.6, N_Cr & BA_tot >0.5, but VIFs look fine
#  Mypairs(Z = MT_ABLAinfo[,comps])#, color=c("red3","purple","blue3")[unclass(CO_ABLAinfo$Band)])
#  corvif(dataz= MT_ABLAinfo[,c("N_Cr","BA_totLR","in5_tot", "ACF")]) 
#  # nothing >0.6 (BA_other & in5_other=0.78, but they're ok), VIFs look fine
#  Mypairs(Z = WA_TSHEinfo[,comps])#, color=c("red3","purple","blue3")[unclass(CO_ABLAinfo$Band)])
#  corvif(dataz= WA_TSHEinfo[,c("N_Cr","BA_totLR","in5_tot", "ACF")]) 
#  # only N_Cr & in5_tot =0.58, but VIFs all look good
#  corvif(dataz= WA_PSMEinfo[,c("N_Cr","BA_totLR","in5_tot", "ACF")]) 
#  Mypairs(Z = WA_PSMEinfo[,comps])#, color=c("red3","purple","blue3")[unclass(CO_ABLAinfo$Band)])
#  # looks good!
#  Mypairs(Z = WA_ABLAinfo[,comps])#, color=c("red3","purple","blue3")[unclass(CO_ABLAinfo$Band)])
#  corvif(dataz= WA_ABLAinfo[,c("N_Cr","BA_totLR","in5_tot", "ACF")]) 
#  # looks fine!

 

######### compselection() function ##############
# now need to find the best model for predicting BAI, so I can get good BAI estimates 
# without the effect of competition or DBH (because both of these vary with elev for most species)
# so I've made a suite of candidate models, and then automated model selection by AIC to get the best
# competition metric for each species. The function takes a Treeinfo dataset and spits out
# the delta AIC table for all models <2 AIC from best model. Then I'm using the ensemble to figure out 
# which variables are generally the best.
# it also spits out the summary of the best model
# and returns the model object for the best model

# things not included in the model selection at present: Band:competition interactions
# perhaps I'll do that after I've narrowed down the field to the best predictors?
# but maybe that fucks me up?
# Also not included: Band:DBH interactions
#      but there only seems to be visual evidence for this relationship to differ in WA & MT TSHE


compselection <- function(dataz, varname="BAI10yr"){
dataz$DBHsc <- scale(dataz$DBH) # scale predictors to avoid numerical instability
dataz$BA_totLRsc <- scale(dataz$BA_totLR) # and to get BAI at species mean predictor level (intercept)
dataz$N_Crsc <- scale(dataz$N_Cr)
dataz$in5_totsc <- scale(dataz$in5_tot)
dataz$ACFsc <- scale(dataz$ACF)
dataz$BAI10yrsc <- scale(dataz$BAI10yr)
null0 <- lmer(get(varname) ~ Band + (1|PairUn), data=dataz, REML=F)
null1 <- lmer(get(varname) ~ Band + DBHsc + (1|PairUn), data=dataz, REML=F)
comp1 <- lmer(get(varname) ~ Band + DBHsc + N_Crsc + (1|PairUn), data=dataz, REML=F)
comp2 <- lmer(get(varname) ~ Band + DBHsc + BA_totLRsc + (1|PairUn), data=dataz, REML=F)
comp3 <- lmer(get(varname) ~ Band + DBHsc + in5_totsc + (1|PairUn), data=dataz, REML=F)
comp4 <- lmer(get(varname) ~ Band + DBHsc + ACFsc + (1|PairUn), data=dataz, REML=F)
comps1 <- lmer(get(varname) ~ Band + DBHsc + BA_sameLR + BA_otherLR + (1|PairUn), data=dataz, REML=F)
comps2 <- lmer(get(varname) ~ Band + DBHsc + in5_same + in5_other + (1|PairUn), data=dataz, REML=F)
compall1 <- lmer(get(varname) ~ Band + DBHsc + N_Crsc + BA_totLRsc + (1|PairUn), data=dataz, REML=F)
compall2 <- lmer(get(varname) ~ Band + DBHsc + N_Crsc + in5_totsc + (1|PairUn), data=dataz, REML=F)
compall3 <- lmer(get(varname) ~ Band + DBHsc + N_Crsc + ACFsc + (1|PairUn), data=dataz, REML=F)
compall4 <- lmer(get(varname) ~ Band + DBHsc + ACFsc + BA_totLRsc + (1|PairUn), data=dataz, REML=F)
compall5 <- lmer(get(varname) ~ Band + DBHsc + ACFsc + in5_totsc + (1|PairUn), data=dataz, REML=F)
compall6 <- lmer(get(varname) ~ Band + DBHsc + BA_totLRsc + in5_totsc + (1|PairUn), data=dataz, REML=F)
compall7 <- lmer(get(varname) ~ Band + DBHsc + BA_totLRsc + in5_totsc + N_Crsc + (1|PairUn), data=dataz, REML=F)
compall8 <- lmer(get(varname) ~ Band + DBHsc + BA_totLRsc + in5_totsc + ACFsc + (1|PairUn), data=dataz, REML=F)
compall9 <- lmer(get(varname) ~ Band + DBHsc + BA_totLRsc + N_Crsc + ACFsc + (1|PairUn), data=dataz, REML=F)
compall10 <- lmer(get(varname) ~ Band + DBHsc + in5_totsc + N_Crsc + ACFsc + (1|PairUn), data=dataz, REML=F)
compall11 <- lmer(get(varname) ~ Band + DBHsc + BA_totLRsc + in5_totsc + N_Crsc + ACFsc + (1|PairUn), data=dataz, REML=F)
compint1 <- lmer(get(varname) ~ Band * N_Crsc + DBHsc + (1|PairUn), data=dataz, REML=F)
compint2 <- lmer(get(varname) ~ Band * BA_totLRsc + DBHsc + (1|PairUn), data=dataz, REML=F)
compint3 <- lmer(get(varname) ~ Band * in5_totsc +  DBHsc + (1|PairUn), data=dataz, REML=F)
compint4 <- lmer(get(varname) ~ Band * ACFsc + DBHsc + (1|PairUn), data=dataz, REML=F)
compint5 <- lmer(get(varname) ~ Band * N_Crsc + BA_totLRsc + Band:BA_totLRsc + DBHsc + (1|PairUn), data=dataz, REML=F)
compint6 <- lmer(get(varname) ~ Band * N_Crsc + in5_totsc + Band:in5_totsc + DBHsc + (1|PairUn), data=dataz, REML=F)
compint7 <- lmer(get(varname) ~ Band * N_Crsc + ACFsc + Band:ACFsc + DBHsc + (1|PairUn), data=dataz, REML=F)
compint8 <- lmer(get(varname) ~ Band * BA_totLRsc + in5_totsc + in5_totsc:Band +  DBHsc + (1|PairUn), data=dataz, REML=F)
compint9 <- lmer(get(varname) ~ Band * BA_totLRsc + ACFsc + ACFsc:Band +  DBHsc + (1|PairUn), data=dataz, REML=F)
compint10 <- lmer(get(varname) ~ Band * in5_totsc +  ACFsc + ACFsc:Band + DBHsc + (1|PairUn), data=dataz, REML=F)

aics <- AIC(null0, null1, comp1, comp2,comp3,comp4,comps1,comps2,compall1, compall2,compall3,compall4,compall5,compall6,compall7,compall8,compall9,compall10,compall11,compint1, compint2, compint3, compint4)
aics$deltaAIC <- aics$AIC - min(aics$AIC)
print(aics[which(aics$deltaAIC <= 2),])
print(summary(get(rownames(aics[which(aics$deltaAIC==0),]))))
return(get(rownames(aics[which(aics$deltaAIC==0),])))
}
#+++++++++++++++++++++++++++++

CO_PIPOinfo$logBAI10yr <- log(CO_PIPOinfo$BAI10yr) # logging goes a good way towards normalizing. but not all the way
CO_PIPOcompmod <- compselection(CO_PIPOinfo, varname="logBAI10yr")
  #shit ton of models that are all similar. ~ACF + in5_tot best model 
  #N_Cr and in5_tot and in5_same+in5_other best single
  # Looks like some toss up between in5_tot, ACF, and N_Cr , with combos all showing up and in5 also showing up with BA_tot
  # looks like two best are ~ACF + in5_tot and ~N_Cr + in5_tot. ACF and N_Cr not highly colinear. in5_tot and N_Cr have cor of 0.69
  # ** Growth INCREASES with Elev
      # once I include one-way interactions, Band * in5_tot is def the best.
CO_PIPOcompCV <- compselection(CO_PIPOinfo, varname="BAI10yrCV")
  # no DBH best, DECREASE with Elev
CO_POTRinfo$logBAI10yr <- log(CO_POTRinfo$BAI10yr) # logging goes a good way towards normalizing. but not all the way
CO_POTRcompmod <- compselection(CO_POTRinfo, varname="logBAI10yr")
  # ~ACF + BA_totLR best model, adding in in5_tot or N_Cr -> similar AICs
  # ** Growth INCREASES with Elev
CO_POTRcompCV <- compselection(CO_POTRinfo, varname="BAI10yrCV")
CO_ABLAinfo$logBAI10yr <- log(CO_ABLAinfo$BAI10yr) # logging goes a good way towards normalizing. but not all the way
CO_ABLAcompmod <- compselection(CO_ABLAinfo, varname="logBAI10yr") # Shit. Got some NAs that I need to fill here.
  # ton of similar models, but best and simplest:
  # ~ BA_totLR
  # ** Growth DECREASES with Elev
CO_ABLAcompCV <- compselection(CO_ABLAinfo, varname="BAI10yrCV")
dataz <- CO_ABLAinfo
test <- update(CO_ABLAcompCV, .~.-Band)
anova(test,CO_ABLAcompCV) # def significant elev effect.

MT_TSHEcompmod <- compselection(MT_TSHEinfo)
test <- update(MT_TSHEcompmod, .~.-Band:N_Crsc)
anova(test, MT_TSHEcompmod)
  # ~N_Cr and anything with N_Cr ** and raw BAI looks similar across elev, but low Elev DBH >> mid and high so actually much slower growing
  # ** Growth INCREASES with elev
MT_TSHEcompCV <- compselection(MT_TSHEinfo, varname="BAI10yrCV")
dataz <- MT_TSHEinfo
test <- update(MT_TSHEcompCV, .~.-Band)
anova(test,MT_TSHEcompCV) # p=0.03041

MT_PSMEinfo$logBAI10yr <- log(MT_PSMEinfo$BAI10yr)
MT_PSMEcompmod <- compselection(MT_PSMEinfo, varname = "logBAI10yr")
  # HUH! BA_same + BA_other best model, but delta AIC over null = 1.46, and LRT over null p=0.06518
  # Only species so far to show a difference between con and heterospecifics, and only just barely.
  #   and not at all in the direction we would think!! (and nothing sinificant marginally in the model itself)
  # ** Growth NS change
dataz <- MT_PSMEinfo # really looks like best to pull out DBH and Band, which means I need to do this the hard way:
ptest1 <- lmer(BAI10yr ~ 1 + (1|PairUn), data=dataz, REML=F)
ptest2 <- lmer(BAI10yr ~ Band + (1|PairUn), data=dataz, REML=F)
ptest3 <- lmer(BAI10yr ~ BA_totLR + (1|PairUn), data=dataz, REML=F)
ptest4 <- lmer(BAI10yr ~ BA_sameLR + BA_otherLR + (1|PairUn), data=dataz, REML=F)
ptest5 <- lmer(BAI10yr ~ BA_sameLR + BA_otherLR + Band + (1|PairUn), data=dataz, REML=F)
ptest6 <- lmer(BAI10yr ~ BA_totLR * Band + (1|PairUn), data=dataz, REML=F)
AIC(ptest1, ptest2,ptest3,ptest4,ptest5, ptest6)
anova(ptest1,ptest4) # p=0.1146, so null model is the best model, no DBH, no BAND
anova(test,MT_PSMEcompmod) # p=0.05229, so marginally significantly better than NULL - DBH
MT_PSMEcompCV <- compselection(MT_PSMEinfo, varname="BAI10yrCV")
step(model = MT_PSMEcompCV,reduce.fixed = T, reduce.random=T)
# so it turns out nothing except N_Cr is significant. best model = CV ~ N_Cr, with 
pCVtest1 <- lmer(BAI10yrCV~ N_Cr + (1|PairUn), data=MT_PSMEinfo, REML=F)
pCVtest2 <- lmer(BAI10yrCV~ BAI10yr + (1|PairUn), data=MT_PSMEinfo, REML=F)
pCVtest3 <- lmer(BAI10yrCV~ BAI10yr + N_Cr + (1|PairUn), data=MT_PSMEinfo, REML=F)
AIC(pCVtest1, pCVtest2, pCVtest3) # adding in BAI10yr makes the model better by 1.83 AIC
anova(pCVtest1, pCVtest3) #LRT p=0.05025


MT_ABLAcompmod <- compselection(MT_ABLAinfo)
dataz <- MT_ABLAinfo
test <- update(MT_ABLAcompmod, .~.-Band)
anova(test,MT_ABLAcompmod) # p=0.08467
  # ~ACF + in5_tot best (both significant), but pretty much anything with ACF showed up as important
  # **Growth DECREASES with Elev, but LRT p only marginally significant
MT_ABLAcompCV <- compselection(MT_ABLAinfo, varname="BAI10yrCV")
  # not perfectly normal residuals (because got some big outliers at the high end)

WA_TSHEinfo$logBAI10yr <- log(WA_TSHEinfo$BAI10yr)
WA_TSHEcompmod <- compselection(WA_TSHEinfo, varname="logBAI10yr")
  # NULL best, but then ACF and N_Cr (and them together) were best of the metrics
  # **Growth NS change
dataz <- WA_TSHEinfo
test <- update(WA_TSHEcompmod, .~.-Band)
anova(test,WA_TSHEcompmod) # p=0.1209 untransformed, p=0.34 transformed
WA_TSHEcompCV <- compselection(WA_TSHEinfo, varname="BAI10yrCV")
dataz <- WA_TSHEinfo
test <- update(WA_TSHEcompCV, .~.-Band)
anova(test,WA_TSHEcompCV) # p=0.2228

WA_PSMEcompmod <- compselection(WA_PSMEinfo)
  # shit ton look good, best: ~BA_tot
  # anything with BA_tot shows up as best model, and BA_same + BA_tot also shows up
  # **Growth INCREASES with Elev
dataz <- WA_PSMEinfo
test <- update(WA_PSMEcompmod, .~.-Band)
anova(test,WA_PSMEcompmod) # p=0.07611 (looks normal)
WA_PSMEcompCV <- compselection(WA_PSMEinfo, varname="BAI10yrCV")
dataz <- WA_PSMEinfo
test <- update(WA_PSMEcompCV, .~.-Band)
test2 <- update(test, .~.+BAI10yr)
test3 <- update(test2, .~.-Band)
anova(test,WA_PSMEcompCV) # p=0.07611 (looks normal)
anova(test, test2) # no help by adding in mean, but not crazy normal...
anova(test2,test3) # elev still marginally significant with BAI10yr added in p=0.06293

WA_ABLAcompmod <- compselection(WA_ABLAinfo)
  # ~ BA_tot + in5_tot
  # anything with in5_tot shows, and then BA_tot + in5_tot + anything shows
  # **Growth greatly INCREASES with Elev
dataz <- WA_ABLAinfo
test <- update(WA_ABLAcompmod, .~.-Band)
anova(test,WA_ABLAcompmod) # p=0.001237 (looks somewhat normal)
WA_ABLAcompCV <- compselection(WA_ABLAinfo, varname="BAI10yrCV")
test <- update(WA_ABLAcompCV, .~.-Band)
test2 <- update(WA_ABLAcompCV, .~.+BAI10yr)
anova(test, WA_ABLAcompCV) # p=0.01589
anova(WA_ABLAcompCV, test2) # damn, adding in BAI10yr makes for a better model (p=0.0184)
plot(fitted(WA_ABLAcompCV)~BAI10yr, WA_ABLAinfo)
plot(resid(WA_ABLAcompCV)~BAI10yr, WA_ABLAinfo)
qqnorm(resid(WA_ABLAcompCV));qqline(resid(WA_ABLAcompCV)) # double damn. Very much nonnormal. Very much shows a relationship with BAI10yr
qqnorm(resid(test2)); qqline(resid(test2)) # but STILL not super normal with BAI10yr added in.

###### Stacked CV, mean BAI for 1 mountain

#### COLORADO _________________________________
quartz(width=6, height=4)
par(mfcol=c(2,3), mar=c(0,3,0,0), oma=c(3,2,2,1))
xvals <- c(0.6666,1.66666,2.6666)
colors <- brewer.pal(n=5, "Set1")[c(4,5,3)] 
# old CV ylim calcs: ylim=c(min(CVgr[,1])-0.3*min(CVgr[,1]),max(CVgr[,1])+0.3*max(CVgr[,1])),

### CO_PIPO ##
meanmod <- CO_PIPOcompmod
CVmod <- CO_PIPOcompCV
dataz <- CO_PIPOinfo
meanmod1 <- update(meanmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
CVgr$xval <- xvals
plot(CVgr[,1]~xvals, xlim=c(0,3.333), bty="n", xaxt="n", ylab="CV in growth", col=colors[1], ylim=c(0.15,0.4)) 
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = CVgr,length = 0, lwd=2, col=colors[1]) 
mtext(text = "CO_PIPO",side = 3)
mtext(text="Growth CV", side=2, line=2)
barplot(height = meangr[,1], xlim=c(0,3.33), space=0.5, width = 0.66666, ylim=c(0, max(meangr[,1])+0.2*max(meangr[,1])), ylab="Mean Annual Growth (mm2)", col=paste0(colors[1],"66"), border = colors[1])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = meangr,length = 0, lwd=2, col=colors[1])
axis(1,at=xvals, labels = c("L","M","H"))
mtext(text=expression(paste("Mean Annual Growth (",mm^2,")", sep="")),side=2, line=2)



### CO_POTR ##
meanmod <- CO_POTRcompmod
CVmod <- CO_POTRcompCV
dataz <- CO_POTRinfo
meanmod1 <- update(CO_POTRcompmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
CVgr$xval <- xvals
plot(CVgr[,1]~xvals, xlim=c(0,3.333), ylim=c(0.15,0.4), bty="n", xaxt="n", ylab="CV in growth", col=colors[2])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = CVgr,length = 0, lwd=2, col=colors[2]) 
mtext(text = "CO_POTR",side = 3)
barplot(height = meangr[,1], xlim=c(0,3.33), space=0.5, width = 0.66666, ylim=c(0, max(meangr[,1])+0.2*max(meangr[,1])), ylab="Mean Annual Growth (mm2)", col=paste0(colors[2],"66"), border = colors[2])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = meangr,length = 0, lwd=2, col=colors[2])
axis(1,at=xvals, labels = c("L","M","H"))

### CO_ABLA ##
meanmod <- CO_ABLAcompmod
CVmod <- CO_ABLAcompCV
dataz <- CO_ABLAinfo
meanmod1 <- update(meanmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
CVgr$xval <- xvals
plot(CVgr[,1]~xvals, xlim=c(0,3.333), ylim=c(0.15,0.4), bty="n", xaxt="n", ylab="CV in growth", col=colors[3])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = CVgr,length = 0, lwd=2, col=colors[3]) 
mtext(text = "CO_ABLA",side = 3)
barplot(height = meangr[,1], xlim=c(0,3.33), space=0.5, width = 0.66666, ylim=c(0, max(meangr[,1])+0.2*max(meangr[,1])), ylab="Mean Annual Growth (mm2)", col=paste0(colors[3],"66"), border = colors[3])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = meangr,length = 0, lwd=2, col=colors[3])
axis(1,at=xvals, labels = c("L","M","H"))



#### MONTANA _________________________________
quartz(width=6, height=4)
par(mfcol=c(2,3), mar=c(0,3,0,0), oma=c(3,2,2,1))
xvals <- c(0.6666,1.66666,2.6666)
colors <- brewer.pal(n=3, "Set1") 
# old CV ylim calcs: ylim=c(min(CVgr[,1])-0.3*min(CVgr[,1]),max(CVgr[,1])+0.3*max(CVgr[,1])),

### MT-TSHE ##
meanmod <- MT_TSHEcompmod
CVmod <- MT_TSHEcompCV
dataz <- MT_TSHEinfo
meanmod1 <- update(meanmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
CVgr$xval <- xvals
plot(CVgr[,1]~xvals, xlim=c(0,3.333), bty="n", xaxt="n", ylab="CV in growth", col=colors[1], ylim=c(0.15,0.4)) 
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = CVgr,length = 0, lwd=2, col=colors[1]) 
mtext(text = "MT-TSHE",side = 3)
mtext(text="Growth CV", side=2, line=2)
barplot(height = meangr[,1], xlim=c(0,3.33), space=0.5, width = 0.66666, ylim=c(0, max(meangr[,1])+0.2*max(meangr[,1])), ylab="Mean Annual Growth (mm2)", col=paste0(colors[1],"66"), border = colors[1])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = meangr,length = 0, lwd=2, col=colors[1])
axis(1,at=xvals, labels = c("L","M","H"))
mtext(text=expression(paste("Mean Annual Growth (",mm^2,")", sep="")),side=2, line=2)

### MT-PSME ##
meanmod <- MT_PSMEcompmod
CVmod <- MT_PSMEcompCV
dataz <- MT_PSMEinfo
meanmod1 <- update(CO_POTRcompmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
CVgr$xval <- xvals
plot(CVgr[,1]~xvals, xlim=c(0,3.333), ylim=c(0.15,0.4), bty="n", xaxt="n", ylab="CV in growth", col=colors[2])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = CVgr,length = 0, lwd=2, col=colors[2]) 
mtext(text = "MT_PSME",side = 3)
barplot(height = meangr[,1], xlim=c(0,3.33), space=0.5, width = 0.66666, ylim=c(0, max(meangr[,1])+0.2*max(meangr[,1])), ylab="Mean Annual Growth (mm2)", col=paste0(colors[2],"66"), border = colors[2])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = meangr,length = 0, lwd=2, col=colors[2])
axis(1,at=xvals, labels = c("L","M","H"))

### MT_ABLA ##
meanmod <- MT_ABLAcompmod
CVmod <- MT_ABLAcompCV
dataz <- MT_ABLAinfo
meanmod1 <- update(meanmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
CVgr$xval <- xvals
plot(CVgr[,1]~xvals, xlim=c(0,3.333), ylim=c(0.15,0.4), bty="n", xaxt="n", ylab="CV in growth", col=colors[3])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = CVgr,length = 0, lwd=2, col=colors[3]) 
mtext(text = "MT_ABLA",side = 3)
barplot(height = meangr[,1], xlim=c(0,3.33), space=0.5, width = 0.66666, ylim=c(0, max(meangr[,1])+0.2*max(meangr[,1])), ylab="Mean Annual Growth (mm2)", col=paste0(colors[3],"66"), border = colors[3])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = meangr,length = 0, lwd=2, col=colors[3])
axis(1,at=xvals, labels = c("L","M","H"))




#### WASHINGTON _________________________________
quartz(width=6, height=4)
par(mfcol=c(2,3), mar=c(0,3,0,0), oma=c(3,2,2,1))
xvals <- c(0.6666,1.66666,2.6666)
colors <- brewer.pal(n=3, "Set1") 
# old CV ylim calcs: ylim=c(min(CVgr[,1])-0.3*min(CVgr[,1]),max(CVgr[,1])+0.3*max(CVgr[,1])),

### WA-TSHE ##
meanmod <- WA_TSHEcompmod
CVmod <- WA_TSHEcompCV
dataz <- WA_TSHEinfo
meanmod1 <- update(meanmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:2,1:2])
meangr$xval <- c(1.66666,2.66666)
CVgr <- data.frame(summary(CVmod1)$coefficients[1:2,1:2])
CVgr$xval <- c(1.66666,2.66666)
plot(CVgr[,1]~c(1.66666,2.66666), xlim=c(0,3.333), bty="n", xaxt="n", ylab="CV in growth", col=colors[1], ylim=c(0.15,0.4)) 
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = CVgr,length = 0, lwd=2, col=colors[1]) 
mtext(text = "WA_TSHE",side = 3)
mtext(text="Growth CV", side=2, line=2)
barplot(height = c(0,meangr[,1]), xlim=c(0,3.33), space=0.5, width = 0.66666, ylim=c(0, max(meangr[,1])+0.2*max(meangr[,1])), ylab="Mean Annual Growth (mm2)", col=paste0(colors[1],"66"), border = colors[1])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = meangr,length = 0, lwd=2, col=colors[1])
axis(1,at=xvals, labels = c("L","M","H"))
mtext(text=expression(paste("Mean Annual Growth (",mm^2,")", sep="")),side=2, line=2)

### MT-PSME ##
meanmod <- WA_PSMEcompmod
CVmod <- WA_PSMEcompCV
dataz <- WA_PSMEinfo
meanmod1 <- update(CO_POTRcompmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
CVgr$xval <- xvals
plot(CVgr[,1]~xvals, xlim=c(0,3.333), ylim=c(0.15,0.4), bty="n", xaxt="n", ylab="CV in growth", col=colors[2])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = CVgr,length = 0, lwd=2, col=colors[2]) 
mtext(text = "WA_PSME",side = 3)
barplot(height = meangr[,1], xlim=c(0,3.33), space=0.5, width = 0.66666, ylim=c(0, max(meangr[,1])+0.2*max(meangr[,1])), ylab="Mean Annual Growth (mm2)", col=paste0(colors[2],"66"), border = colors[2])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = meangr,length = 0, lwd=2, col=colors[2])
axis(1,at=xvals, labels = c("L","M","H"))

### WA_ABLA ##
meanmod <- WA_ABLAcompmod
CVmod <- WA_ABLAcompCV
dataz <- WA_ABLAinfo
meanmod1 <- update(meanmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
CVgr$xval <- xvals
plot(CVgr[,1]~xvals, xlim=c(0,3.333), ylim=c(0.15,0.4), bty="n", xaxt="n", ylab="CV in growth", col=colors[3])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = CVgr,length = 0, lwd=2, col=colors[3]) 
mtext(text = "WA_ABLA",side = 3)
barplot(height = meangr[,1], xlim=c(0,3.33), space=0.5, width = 0.66666, ylim=c(0, max(meangr[,1])+0.2*max(meangr[,1])), ylab="Mean Annual Growth (mm2)", col=paste0(colors[3],"66"), border = colors[3])
error_bars(xvar ="xval", yvar="Estimate",upper="Std..Error",errordata = meangr,length = 0, lwd=2, col=colors[3])
axis(1,at=xvals, labels = c("L","M","H"))







############ Competition Analysis #################
# now to actually extract the importance of competition,
# and see if it changes across elevation in a predictable way. 



plot(BAI10yr~DBH, WA_TSHEinfo, pch=16, col=Band)
