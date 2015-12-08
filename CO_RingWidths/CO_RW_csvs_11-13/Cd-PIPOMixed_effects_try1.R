#=======================================================
#-------------------------------------------------------
#               RW mixed effects modeling
#               Leander DL Anderegg
#                   11/11/13
#-------------------------------------------------------
#=======================================================

## load core.average() function


rm(list=ls())
setwd("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13")
# set working directory to my scratch file in Dropbox
library(dplR)
require(reshape)
require(dplyr)

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


##################################################
#       read in core and climate data 
#       (if you don't need to change the process)
##################################################

# with raw climate values
rw.clim.allraw <- read.csv("CO-PIPO-allRWI_rawClim.csv")

# with scaled climate values
rw.clim.allelev <- read.csv("CO-PIPO-allRWIClim.csv")

# PIPO master chronologies
PIPO.master <- read.csv("MASTER_CO-PIPO_02072014.csv")




#==============================================================================
#===================== Code used to create the above .csvs ====================

################################
#     read tree core data
################################
## low elevation
csv1 <- "CO-PIPO-L_RWdata11_7_13test2.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
treedata1 <- core_average(data=data1)
## mid elevation
csv2 <- "CO-PIPO-M_RWdata11_7_13.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "CO-PIPO-H_RWdata10_31_13.csv"
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3)


### Note: in this formulation, treedata [,5:115] = yrs 2011-1902
# theoretically making everything 2011-1902, 
# climate goes 1901-2011, but Low elevation only goes to 1902
#treedata <- rbind(treedata1[,1:115], treedata2[,1:115], treedata3[,1:115])


###########################################
#           Clean tree core Data 
###########################################

#========= subsetting to desired dates============
#-------------------------------------------------

#______________________________________
### customize these to get dates
# for PRISM data from climate WNA, first full water year = 1902, last year =2011
#workingdata <-treedata1
start <- "1902"   # first date for time series
end <- "2011"     # last date of time series
method <- "Spline" # detrending method desired. c("Spline", "ModNegExp", "Mean")
dropyrs <- 3
####

#### Making both a dataframe of detrended rwi and trend-in ring widths
#rw.det <- core_trunc(workingdata=workingdata, detmethod=method, detrend=TRUE)
#rw.nodet <- core_trunc(workingdata=workingdata, detmethod=method, detrend=FALSE)

# for specific elevations
# High Elevation
rw.detH <- core_trunc(workingdata=treedata3, first=start, last=end, detmethod=method, detrend=TRUE, drop=dropyrs)
# Mid Elevation
rw.detM <- core_trunc(workingdata=treedata2, first=start, last=end,  detmethod=method, detrend=TRUE, drop=dropyrs)
# Low Elevation
rw.detL <- core_trunc(workingdata=treedata1, first=start, last=end,  detmethod=method, detrend=TRUE, drop=dropyrs)


### Making a .csv of master chronologies
  # this just uses the mean.
L.master <- apply(X=rw.detL, MARGIN=1,FUN=mean, na.rm=TRUE)
M.master <- apply(X=rw.detM, MARGIN=1,FUN=mean, na.rm=TRUE)
H.master <- apply(X=rw.detH, MARGIN=1,FUN=mean, na.rm=TRUE)
  # or using the biweighted robust mean (more the industry standard)
L.masterbw <- chron(rw.detL, prefix="L")
M.masterbw <- chron(rw.detM, prefix="M")
H.masterbw <- chron(rw.detH, prefix="H")

PIPO.master <- data.frame(L.master, M.master, H.master)
names(PIPO.master) <- c("PIPO-L", "PIPO-M", "PIPO-H")
#write.csv(PIPO.master, file="MASTER_CO-PIPO_02072014.csv")

#--------------------------------------
#==============================================




#########################################
#        Read and process climate data
#######################################

#==============================================
#------------- Climate data for analysis ------

#-----------  PIPO
# clim <- read.csv("PIPO_climate.csv") # this is just monthly values
clim <- read.csv("PIPO_climate_Annetc.csv") 

#----------- HIGH ELEV
climH <- clim[which(clim$ID2=="PIPO-H"),-c(2:6)]
climH <- climH[-1,] # get rid of 1901 because not a complet water year
library(plyr)
climH.des <- arrange(climH, desc(Year))
# make sure everything's numeric
for (i in 2:ncol(climH.des)){ climH.des[,i] <- as.numeric(as.character(climH.des[,i]))} 
Year <- climH.des$Year

# scale climate variables to limit covariance in random effects model
climH.des.scale <- data.frame(Year, apply(X=climH.des[,-1],MARGIN=2,FUN=scale))

#----------Mid Elevation
climM <- clim[which(clim$ID2=="PIPO-H"),-c(2:6)]
climM <- climM[-1,] # get rid of 1901 because not a complet water year
climM.des <- arrange(climM, desc(Year))
# make sure everything's numeric
for (i in 2:ncol(climM.des)){ climM.des[,i] <- as.numeric(as.character(climM.des[,i]))} 
Year <- climM.des$Year

# scale climate variables to limit covariance in random effects model
climM.des.scale <- data.frame(Year, apply(X=climM.des[,-1],MARGIN=2,FUN=scale))

#----------Low Elevation
climL <- clim[which(clim$ID2=="PIPO-H"),-c(2:6)]
climL <- climL[-1,] # get rid of 1901 because not a complet water year
climL.des <- arrange(climL, desc(Year))
# make sure everything's numeric
for (i in 2:ncol(climL.des)){ climL.des[,i] <- as.numeric(as.character(climL.des[,i]))} 
Year <- climL.des$Year

# scale climate variables to limit covariance in random effects model
climL.des.scale <- data.frame(Year, apply(X=climL.des[,-1],MARGIN=2,FUN=scale))


############# Scaled Climate Data Frames #################
climH.des.scale <- data.frame(Year, apply(X=climH.des[,-1],MARGIN=2,FUN=scale))
climM.des.scale <- data.frame(Year, apply(X=climM.des[,-1],MARGIN=2,FUN=scale))
climL.des.scale <- data.frame(Year, apply(X=climL.des[,-1],MARGIN=2,FUN=scale))
#########################



#__________________________________________________________________
# Stuff for trying to determine best climate variables ------------

###### VARIANCE INFLATION FACTOR TEST ####################
# Trying to assess varaince inflation factors from Zuur "Mixed Effects Models"
  # but currently doesn't work if I have more than 7 variables in it
source("/Users/leeanderegg/Desktop/ZuurMixedModelling/AllRCode/HighstatLibV6.R")
corvif(climH.des[,c("TaveGS", "PPTGS", "PPTAn", "DD5_An","NFFDGS", "CMDGS", "CMDAn")])
## This suggests that only PPTAn and NFFDGS are <6, with DD5_An and CMDAn <15

#######
climH.test <- as.data.frame(climH.des[,-1])
pairs(climH.test, lower.panel=panel.smooth2,
      upper.panel=panel.cor,diag.panel=panel.hist)
####### hmm. many things are strongly colinear. don't know the best way to remove them

# temperature related variables:
tempvars <- c("TmaxAn", 'TminAn', "TaveAn", "TaveWin", "TaveGS", "DD5_An", "NFFDGS", "NFFDAn", "ErefGS")
quartz(width=5, height=5)
pairs(climH.des[,tempvars], lower.panel=panel.smooth2,
      upper.panel=panel.cor,diag.panel=panel.hist)
  # from this plot I conclude that TaveWin, TaveGS, and TmaxAn + TminAn
    # are probably the ones to stick with. TaveGS and all of the derived variables
    # are extremely correlated. So I should only include one of either
    # TmaxAn, TaveGS, or ErefGS (but if ErefGS, don't use CMDGS, use a different precip)
corvif(climH.des[,c("TaveGS", "DD5_An","NFFDGS", "NFFDAn", "ErefGS")])
  # huh. don't know how to interpret this.
## Possible T combos: ##
  # TmaxAn + TminAn, (TmaxAn * TminAn)
  # TaveGS + TaveWin, (TaveGS * TaveWin)
  # TaveGS + TminAn, TmaxAn + TaveWin
  # ErefGS + TminAn
  # ErefGS + TaveWin
  # TaveAn

# precipitation related variables:
pptvars <- c("PPTGS", "PPTAn", "PPTWin", "PASAn", "CMDGS", "CMDAn", "TaveAn")
quartz(width=5, height=5)
pairs(climH.des[,pptvars], lower.panel=panel.smooth2,
      upper.panel=panel.cor,diag.panel=panel.hist)
corvif(climH.des[,pptvars]) ## doesn't work and I don't know why,
                             # seems to never have tmp_cor in myvif called by corvif
  # seems like PPTGS is correlated with and driving both CMDs, so probalby just use PPTGS?
  # don't use PPTAn, use one summer PPT (PPTGS, CMDGS, or CMDAn) and one winter (PPTWin/PASAn)
    # also, TaveGS & TmaxAn are not strongly correlated with CMDGS or CMDAn, though ErefGS is
## Possible PPT combos: ##
  # PPTGS + PPTWin
  # PPTGS + PASAn
  # CMDGS + PPTWin
  # CMDGS + PASAn
  # CMDAn + PPTWin
  # CMDAn + PASAn
  # PPTAn
  # Relevant interactions: Summer T and Summer Precip
#_____________________________________________________________________________________


#==========================================================
###############  Put together full rwi and climate data
#=====================================================
# Put together the full data frame with rw and clim
combine_rwclim <- function(rw.det,clim.des.scale,elev,Treename="Tree") {
  rw.clim <- data.frame(rw.det,clim.des.scale)
  id.vars <- colnames(clim.des.scale)
  measure.vars <- colnames(rw.det) # assuming you want detrended rwi

require(reshape)
  rw.clim.df <- melt (rw.clim, id.vars = id.vars,variable_name=Treename)
  rw.clim.df$Elev <- as.factor(rep(elev, times=nrow(rw.clim.df)))

# adding in competition factor
### Breaking apart H and L
  # note, have to change whether grepping on H or L based on elevation (grep on H for L elev)
  if(elev=="H"){
    rw.clim.N <- rw.clim.df[grep(pattern="L",x=rw.clim.df$Tree),] 
    rw.clim.C <- rw.clim.df[-grep(pattern="L",x=rw.clim.df$Tree),]
    rw.clim.N$Comp <- as.factor(rep("N", times = nrow(rw.clim.N)))
    rw.clim.C$Comp <- as.factor(rep("C", times= nrow(rw.clim.C)))
  }
  else {
    rw.clim.C <- rw.clim.df[grep(pattern="H",x=rw.clim.df$Tree),] 
    # note, have to change whether grepping on H or L based on elevation (grep on H for L elev)
    rw.clim.N <- rw.clim.df[-grep(pattern="H",x=rw.clim.df$Tree),]
    rw.clim.N$Comp <- as.factor(rep("N", times = nrow(rw.clim.N)))
    rw.clim.C$Comp <- as.factor(rep("C", times= nrow(rw.clim.C)))
  }
#       ======= data frames to use ==========    #
#rw.clim.df # detrended rw and climate data by tree and year and ready to go
  rw.clim.df.narm <- na.omit(rbind(rw.clim.C, rw.clim.N)) # same df but na's removed 'casue {nlme} doesn't like them
return(rw.clim.df.narm)
}
#=======================================================
# for dealing with scaled climate variables
Hrw.clim.df.narm <-combine_rwclim(rw.det=rw.detH,clim.des.scale=climH.des.scale,elev="H")
Mrw.clim.df.narm <-combine_rwclim(rw.det=rw.detM,clim.des.scale=climM.des.scale,elev="M")
Lrw.clim.df.narm <-combine_rwclim(rw.det=rw.detL,clim.des.scale=climL.des.scale,elev="L")


rw.clim.allelev <- rbind(Hrw.clim.df.narm, Mrw.clim.df.narm, Lrw.clim.df.narm)
write.csv(rw.clim.allelev, file="CO-PIPO-allRWIClim.csv")


# or for dealing with raw climate variables
Hrw.clim.raw <-combine_rwclim(rw.det=rw.detH,clim.des.scale=climH.des,elev="H")
Mrw.clim.raw <-combine_rwclim(rw.det=rw.detM,clim.des.scale=climM.des,elev="M")
Lrw.clim.raw <-combine_rwclim(rw.det=rw.detL,clim.des.scale=climL.des,elev="L")

rw.clim.allraw <- rbind(Hrw.clim.raw, Mrw.clim.raw, Lrw.clim.raw)
#write.csv(rw.clim.allraw, file="CO-PIPO-allRWI_rawClim.csv")



#==================================================
##########     Making quick spark plots of data
#==================================================
#### Change this to get different elevations
toplot <- subset(x=rw.clim.allelev,subset=Elev=="L")
toplot <- toplot[order(toplot$Tree),]
# with H & L layered #THIS IS NOT WORKING, AND I DON'T KNOW WHY
sparkplots(years=toplot$Year,names=toplot$Tree, data=toplot$PPTGS, types=toplot$Comp,legend=TRUE)
# with all plotted as lines
sparkplots(years=toplot$Year,names=toplot$Tree, data=toplot$value, types=rep(x=1,times=nrow(toplot)),legend=FALSE)



#=========================== End Code for date frame creation ====================

#=====================================================
###### Checking for temporal autocorrelation ---------
#=====================================================

#currently just doing this for master chronologies. need to do this for each tree?
# PIPO master chronologies
PIPO.master <- read.csv("MASTER_CO-PIPO_02072014.csv")
acf(PIPO.master$PIPO.L) # 1 and 2 year lag just barely significant
acf(PIPO.master$PIPO.M) # 1 year lag pretty significant
acf(PIPO.master$PIPO.H) # 1 year lag somewhat significant

for(i in unique(rw.clim.df.narm$Tree)){
  tmp <- subset(x=rw.clim.df.narm$value,subset=rw.clim.df.narm$Tree==i)
  quartz(width=5, height=4)
  acf(tmp, main=i)
}

#======================================================
########## RANDOM EFFECTS SELECTION ###################
#======================================================
# Zuur et al. "Mixed Effects Models" suggests a process of
  # 1) Select random effects structures using REML and a "beyond" optimal
      # fixed effects model so that you don't accidentally get fixed effects 
      # stuff in the random effects structure. They use AIC to select
  # 2) Select fixed structure using backwards selection and ML
  # 3) rerun optimal model with REML and validate

# read in data frame of all PIPO detrended RWI and associated climate data
rw.clim.allelev <- read.csv("CO-PIPO-allRWIClim.csv", header=T)
rw.clim.df.narm <- subset(rw.clim.allelev,subset=Elev=="L") # subset to desired elevation

### Random Effects Formulation
  # need to update when I winnow down the clim variables,
  # but for now will just use my judgement
library(lattice)

library("nlme") #this doesn't work with NAs. gives p value (but suspect)
  # also doesn't seem to work with multiple random effects
test1 <- lme(value~PPTAn, data=rw.clim.df.narm, random=~1|Tree)

library("lme4") # this works with NAs. no p value
test2 <- lmer(value~PPTAn + (1|Tree), data=rw.clim.df.narm)

#### formulation for full model: TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn
      # this should cover one variable from almost all relavent categories



##############################################################
#             Random Effects Model Selection
##############################################################

#=== no random effects
r1 <- gls(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp, data = rw.clim.df.narm) 
#=== tree random intercept
#r2 <- lme(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn, random = ~1|Tree, data = rw.clim.df.narm)
r2r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp+(1|Tree), data = rw.clim.df.narm)
#=== tree random slope & intercept, can't get lme to work with necessary random effect structure
#r3 <- lme(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn, random = ~1 + CMDGS|Tree, data = rw.clim.df.narm)
#r3.1 <- lme(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn, random = ~1 + TmaxAn|Tree, data = rw.clim.df.narm)
#r3r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn + (1 + CMDGS|Tree), data = rw.clim.df.narm)
r3r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp 
             +(CMDGS|Tree)+(TmaxAn|Tree)+(TminAn|Tree)+(TaveWin|Tree)+(NFFDGS|Tree) +(ErefGS|Tree)+(PASAn|Tree)
             , data = rw.clim.df.narm)
r3r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp 
            +(ErefGS|Tree)+(Comp|Tree)
            , data = rw.clim.df.narm)
#=== tree random slope & intercept + year random intercept
#r4 <- lme(value~TmaxAn+TminAn+TaveWin+PPTGS+PPTWin+NFFDGS+PASAn, random = ~1|Year+PPTGS|Tree, data = rw.clim.df.narm)
#r4r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn + (CMDGS|Tree) + (1|Year), data = rw.clim.df.narm)
r4r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp 
                    + (0+CMDGS|Tree)+(0+TmaxAn|Tree)+(0+TminAn|Tree)+(0+TaveWin|Tree)+(0+NFFDGS|Tree)+(0+ErefGS|Tree)+(0+PASAn|Tree)+(0+Comp|Tree)
                    + (1|Year)
                    , data = rw.clim.df.narm, control=lmerControl(optCtrl=list(maxfun=10000) ))
r4r1 <- lmer(value~TmaxAn+NFFDGS+ErefGS+CMDGS+PASAn+Comp 
             +(0+TmaxAn|Tree)+(0+NFFDGS|Tree)+(0+ErefGS|Tree)+(0+CMDGS|Tree)+(0+PASAn|Tree)
            + (1|Year)
            , data = rw.clim.df.narm, control=lmerControl(optCtrl=list(maxfun=10000) ))
# tree random slope & intercept + year random slope + intercept
r5r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp 
            + (TmaxAn|Tree)+(TminAn|Tree)+(TaveWin|Tree)+(NFFDGS|Tree)+(ErefGS|Tree)+(CMDGS|Tree)+(PASAn|Tree)+(Comp|Tree)
            + (TmaxAn|Year)+(TminAn|Year)+(TaveWin|Year)+(NFFDGS|Year)+(ErefGS|Year)+(CMDGS|Year)+(PASAn|Year)+(Comp|Year)
            , data = rw.clim.df.narm)
# tree random intercept + year random slope + intercept
r6r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp 
            + (1|Tree) 
            + (TmaxAn|Year)+(TminAn|Year)+(TaveWin|Year)+(NFFDGS|Year)+(ErefGS|Year)+(CMDGS|Year)+(PASAn|Year)+(Comp|Year)
            , data = rw.clim.df.narm)
# tree random intercept + year random intercept
r7r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp + (1|Tree) + (1|Year), data = rw.clim.df.narm)
# year random slope + random intercept
r8r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp 
            + (TmaxAn|Year)+(TminAn|Year)+(TaveWin|Year)+(NFFDGS|Year)+(ErefGS|Year)+(CMDGS|Year)+(PASAn|Year)+(Comp|Year)
            , data = rw.clim.df.narm)
# year random intercept
r9r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp + (1|Year), data = rw.clim.df.narm)
r9 <- lme(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp, random= ~ 1|Year, data = rw.clim.df.narm)
AIC (r1, r2r, r3r, r4r, r5r, r7r, r8r, r9r)
AIC(r1, r2r, r3r, r4r, r7r, r9r) # just the models that don't have Year random slopes
    # Year random intercpt has the best AIC (AIC=-1052.5071) without comp (w/ comp AIC=-1042.6615),
    # with tree random intercept as a close second (AIC=-1050.5044) (w/comp AIC=-1040.6589)
    # Tree random slope & intercept and year random slope seemed to be best AIC (-1054.2156)
    # before I realized slope had to be for each fixed variable.
    # adding comp just decreased AIC for all models without changing order

#### For High Elevation #####
#    df        AIC
#r1  10  -517.9388
#r2r 10  -525.5199
#r3r 34  -471.0092
#r4r 35 -1003.1436
#r7r 12 -1040.6589
#r9r 11 -1042.6615

#### For Mid Elevation #####
#   df       AIC
#r1  10  198.0039
#r2r 10  190.8675
#r3r 34  245.7540
#r4r 35 -754.7968
#r7r 12 -785.8046
#r9r 11 -787.8046

#### For Low Elevation #####
#   df       AIC
#r1  10  606.0307
#r2r 10  599.2634
#r3r 34  653.9398
#r4r 35 -410.1393
#r7r 12 -437.6570
#r9r 11 -439.6570


##### if I include all random slopes in the effects (which I think I might need to), then 1|Year is the clear winner

### Checking the residual structure of my "over the top" model with the correct random effects structure
r9r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp + (1|Year), data = rw.clim.df.narm)
resids <- residuals(object=r9r,type="pearson")
fits <- fitted(r9r)
boxplot(resids~Comp, data=rw.clim.df.narm, main="Competition")
plot(resids~Year, data=rw.clim.df.narm)
boxplot(resids~Tree, data=rw.clim.df.narm);abline(h=0, lty=2)
plot(resids~CMDGS, data=rw.clim.df.narm)
plot(resids~TmaxAn, data=rw.clim.df.narm); abline(h=0, lty=2)
plot(resids~fits)

#### And a quick check of nasty temporal autocorrelation
for(i in unique(rw.clim.df.narm$Tree)){
  tmp <- subset(x=resids,subset=rw.clim.df.narm$Tree==i)
  quartz(width=5, height=4)
  acf(tmp, main=i)
}
# it looks like some of these bad boys do have temporal autocorrelation. Makes me think that I shoudl probably try an AR1 model

#----------- Trying to incorporate and AR(1) correlation structure ----------
r9 <- lme(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp, random= ~ 1|Year, data = rw.clim.df.narm)
r9.ar1 <- lme(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp, random= ~ 1|Year, data = rw.clim.df.narm, correlation=corAR1(form=~Year|Tree))



##======================================================
#---------- Checking out Autocorrelation structure---------
#=======================================================
r9r <- lmer(value~TmaxAn+TminAn+TaveWin+NFFDGS+ErefGS+CMDGS+PASAn+Comp + (1|Year), data = rw.clim.df.narm)
resids <- residuals(object=r9r,type="Normalized")





##======================================================
#---------- A weak first attempt at model selection---------
#=======================================================

Null <- lmer(value~1+(1|Year), data=rw.clim.df.narm, REML=F)
water <- lmer(value~CMDGS+(1|Year), data=rw.clim.df.narm, REML = F) # from model r9, it looks like CMDGS is definitely the most significant
wintemp <- lmer(value~TaveWin + (1|Year), data=rw.clim.df.narm, REML=F)
gstemp <- lmer(value~TaveGS + (1|Year), data=rw.clim.df.narm, REML=F)
comp <- lmer(value~Comp + (1|Year), data=rw.clim.df.narm, REML=F)

anova(Null,water,wintemp, gstemp, comp)

# testing out whether a competition interaction is going on
watercomp <- lme(value ~CMDGS*Comp, random=~1|Year, data= rw.clim.df.narm)
    # interaction might be significant (bad p =0.04)
gstcomp <- lme(value ~TaveGS*Comp, random=~1|Year, data= rw.clim.df.narm)
    # interaction definitely not significant




# Selecting best Summer T variable
T1 <- lmer(value ~ PPTGS + PPTWin + TmaxAn + (1|Year), data = rw.clim.df.narm, REML=FALSE)
T2 <- lmer(value ~ PPTGS + PPTWin + TaveGS + (1|Year), data = rw.clim.df.narm, REML=FALSE)
T3 <- lmer(value ~ PPTGS + PPTWin + ErefGS + (1|Year), data = rw.clim.df.narm, REML=FALSE)
AIC(T1, T2, T3) # order = TmaxAn (-850.2576), ErefGS (-847.0353), TaveGS (-841.1842)
T4 <- lmer(value ~ TmaxAn + (1|Year), data = rw.clim.df.narm, REML=FALSE)
T5 <- lmer(value ~ TaveGS + (1|Year), data = rw.clim.df.narm, REML=FALSE)
T6 <- lmer(value ~ ErefGS + (1|Year), data = rw.clim.df.narm, REML=FALSE)
AIC(T4, T5, T6) # order = ErefGS (-839.8555), TmaxAn (-837.1224), TaveGS (-826.7443)
### Hmm. It looks like TaveGS is the clear loser, but ErefGS and TmaxAn are close

#Selecting Winter T variable
T7 <- lmer(value ~ TminAn + (1|Year), data = rw.clim.df.narm, REML=FALSE)
T8 <- lmer(value ~ TaveWin + (1|Year), data = rw.clim.df.narm, REML=FALSE)
AIC(T7, T8) # They're identical. might go for min for physiology's sake

#Selecting Summer Precip variable
P1 <- lmer(value ~ TmaxAn + TminAn + PPTGS + (1|Year), data = rw.clim.df.narm, REML=FALSE)
P2 <- lmer(value ~ TmaxAn + TminAn + CMDGS + (1|Year), data = rw.clim.df.narm, REML=FALSE)
P3 <- lmer(value ~ TmaxAn + TminAn + PPTAn + (1|Year), data = rw.clim.df.narm, REML=FALSE)
P4 <- lmer(value ~ TmaxAn + TminAn + CMDAn + (1|Year), data = rw.clim.df.narm, REML=FALSE)
AIC(P1, P2, P3, P4) # order = CMDGS (-840.9681, or -849.6085 with PPTWin and -842.1381 with PASAn)
                #         PPTGS (-839.9210, or -849.0158 with PPTWin and -841.8935 with PASAn)
                # so CMDGS is slightly better, but typically <1 dAIC
                # and PPTAn is the best of them all by dAIC ~2, and waaay beter than CMDAn

P5 <- lmer(value ~ PPTGS + PPTWin + (1|Year), data = rw.clim.df.narm, REML=FALSE)
P6 <- lmer(value ~ CMDGS + PPTWin + (1|Year), data = rw.clim.df.narm, REML=FALSE)
P7 <- lmer(value ~ PPTGS + (1|Year), data = rw.clim.df.narm, REML=FALSE)
P8 <- lmer(value ~ CMDGS + (1|Year), data = rw.clim.df.narm, REML=FALSE)
P9 <- lmer(value ~ PPTAn + (1|Year), data = rw.clim.df.narm, REML=FALSE)
P10 <- lmer(value ~ CMDAn + (1|Year), data = rw.clim.df.narm, REML=FALSE)
AIC(P5, P6, P7, P8, P9, P10) # order = CMDGS (-845.2607 with PPTWin or -836.9290 without)
                             #         PPTGS (-839.0732 with PPTWin or -830.0728 without)
# PPTWin definitely does better than PASAn with other variables

#### Take home's about Climate Variables:
# Summer T: TmaxAn and ErefGS are close and depend on whether alone or with PPT. TaveGS sucks
# Winter T: TminAn and TaveWin are essentially identical, and TminAn more physiological?
# Summer P: CMDGS is slightly better than PPTGS with other variables, and considerably better alone
# Winter P: PPTWin def better than PASAn
# Annual P: PPTAn acutlally beats CMDAn typically



#======================================================
#         Fixed Effects  Model Selection           
#======================================================

Mn <- lmer(value~1+(1|Tree)+(1|Year), data=rw.clim.df, REML=F)







# ========= Leftovers and Random code

## plotting PIPO for outreach

quartz(width=8, height=4)
par(mar=c(5,5,0,1))
plot(value-.5~Year, data=Hrw.clim.df.narm[which(Hrw.clim.df.narm$Tree=="H.3L"),], type="l", yaxt="n", ylim=c(0,3), lwd=3, col="green3", ylab="")
points(((PPTAn/5)+.5)~Year, data=Hrw.clim.df.narm[which(Hrw.clim.df.narm$Tree=="H.3L"),], type="l", lwd=3,col="red2")
points((value*.9+1)~Year, data=Lrw.clim.df.narm[which(Lrw.clim.df.narm$Tree=="L.6H"),], type="l", lwd=3,col="green")
points(((PPTAn/5)+2)~Year, data=Lrw.clim.df.narm[which(Lrw.clim.df.narm$Tree=="L.6H"),], type="l", lwd=3,col="red2")
mtext(text="Core(green) v ann precip (red) for high (lower) and low (upper) elev",side=3, padj=2)
#points(L.master~as.numeric(names(L.master)), type="l", lwd=3,col="blue")
plot(climH.des.scale$CMDGS ~ climH.des.scale$Year, type="l", lwd=3)

summary(Hrw.clim.df.narm$value[which(Hrw.clim.df.narm$Tree=="H.3L")]) # mean =1

## Making master chronology
L.master <- apply(X=rw.detL, MARGIN=1,FUN=mean, na.rm=TRUE)
plot(climH.des$TaveAn ~ climH.des$Year, type="l", lwd=3)
points(x=2002, y=climH.des$TaveAn[climH.des$Year==2002], pch=15, col="red")
points(L.master~as.numeric(names(L.master)), type="l", lwd=3,col="blue")

climH.des$Year[which(climH.des$PPTAn==min(climH.des$PPTAn))]

quartz(width=7, height=6)
par(mar=c(0,5,0,1), mfrow=c(3,1), oma=c(5,2,1,0))
plot(climH.des$CMDAn ~ climH.des$Year, type="l", lwd=3, ylab="average annual T")
points(x=2002, y=climH.des$CMDAn[climH.des$Year==2002], pch=16, col="red", cex=1.5)
plot(climH.des$PPTAn ~ climH.des$Year, type="l", lwd=3, ylab="anaul precip (water year)")
points(x=2002, y=climH.des$PPTAn[climH.des$Year==2002], pch=16, col="red", cex=1.5)
plot(climH.des$PPTGS ~ climH.des$Year, type="l", lwd=3, ylab="Apr-Sept Precip")
points(x=2002, y=climH.des$PPTGS[climH.des$Year==2002], pch=16, col="red", cex=1.5)



Null <- lmer(value~1+(1|Tree)+(1|Year), data=rw.clim.df, REML=F)
water <- lmer(value~CMDGS+(1|Tree)+(1|Year), data=rw.clim.df, REML = F)
temp <- lmer(value~TmaxAn + (1|Tree) + (1|Year), data=rw.clim.df, REML=F)
anova(Null,water,temp)

quartz(width=10, height=6)
par(mfrow=c(3,5))
for (i in 2:16) plot(rw.clim.df$value~rw.clim.df[,i], ylab = "RWI", xlab = colnames(rw.clim.df)[i])

quartz(width=5, height=6)
plot(value~PPTAn, data = rw.clim.df, col = Tree)
abline(lm(value~PPTAn, data=rw.clim.df), lwd=3)
trees <- unique(rw.clim.df$Tree)
for (i in trees) abline(lm(value~1+PPTAn, data=rw.clim.df[rw.clim.df$Tree==i,]), lty=2)


test <- lmer(value~CMDGS + (CMDGS|Tree) + (1|Year), data=rw.clim.df.narm)
effects <- ranef(test)
effects$Tree[,2]
coef(test)

plot(value~TmaxAn, data=rw.clim.df)



#-----------------------------
#-------------------------------------------
#------------------------------------------------------------------











treedata <-core_average(data=data3)
plot(y=t(treedata[1,5:112]),x=colnames(treedata)[5:112], type="l")
lines(y=t(treedata[2,5:112]),x=colnames(treedata)[5:112])
lines(y=t(treedata[5,5:112]),x=colnames(treedata)[5:112])
lines(y=t(treedata[8,5:112]),x=colnames(treedata)[5:112])
points(y=treedata[5,"1986"], x = 1986)
points(y=treedata[1,"1987"], x = 1987)
points(y=treedata[1,"1988"], x= 1988)
### now detrending data




