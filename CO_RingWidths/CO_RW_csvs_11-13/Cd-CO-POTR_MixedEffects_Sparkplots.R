#=======================================================
#-------------------------------------------------------
#               RW mixed effects modeling (for POTR)
#               Leander DL Anderegg
#                   11/11/13
#-------------------------------------------------------
#=======================================================

## load core.average() function


rm(list=ls())
setwd("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13")
setwd("C:/Users/Western Hemlock/Dropbox/Postdocs&GradProjects/TreeCoresRanges(Leander)/CO_RingWidths/CO_RW_csvs_11-13/")
# set working directory to my scratch file in Dropbox
library(dplR)

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
## For Western Hemclock
#source("C:/Users/Western Hemlock/Dropbox/Postdocs&GradProjects/TreeCoresRanges(Leander)/CO_RingWidths/CO_RW_csvs_11-13/fn_core_trunc.R")
#source("C:/Users/Western Hemlock/Dropbox/Postdocs&GradProjects/TreeCoresRanges(Leander)/CO_RingWidths/CO_RW_csvs_11-13/fn_core_average.R")

#############################
#     Function to keep each core unaveraged
#############################
core_clean <- function(data, yrcored = 2013)
{
  data$Analysis.Date.Time <- as.POSIXct(strptime(data$Analysis.Date.Time, format = "%m/%d/%Y %H:%M"))
  un <- unique(data$TreeID)
  dataclean <- data[1:length(un),]
  for(i in 1:length(un))
  {
    data_subset <- data[which(data$TreeID==un[i]),]
    dataclean[i,] <- data_subset[which(data_subset$Analysis.Date.Time  == max(data_subset$Analysis.Date.Time)),]
  }
  #dataclean$Tag
  
  coredataraw <- dataclean
  
  
  badcols <- c( "Path.identification","Site.identification" ,"YearLastRing", "Sapwood","Tree.height","Tree.age"
                ,"SectionHeight" , "User.variable" ,"DataType","OffsetToNext","ImageName"
                , "Acquisition.Date.Time","Modified.Date.Time", "ImageSizeH.V.NBits.Channel"
                ,"CalibMethod.XCal.YCal.EditedDendro","ScannerCamera.Make.Model.Software","LensFocLength..35mm."
                ,"PathBegX.BegY.EndX.EndY.Width","RingBoundary.AutoMan.Meth.Precise","EarlywoodDef"
                ,"DensActive.Media.Calib","DensNSteps.MatDens.Interpol.SampleThick","DensStepsThick"
                ,"DensStepsLightInt","DensStepsWoodDens","DiskArea","DiskPerim","DiskAvgDiam"
                ,"DiskFormCoef","CompWoodArea","VoidArea","PathLength")
  # columns to kill from Windendro auto output
  coredata <- coredataraw[,!(colnames(coredataraw) %in% badcols)]
  #coredata[,"X"] <- NULL # removing 2013 (year core was collected)
  #  return(coredata)
  trees <- unique(coredata$TreeID)
  Transect <- rep(coredata$Transect[1], times=length(trees))
  Species <- rep(coredata$Species[1], times=length(trees))
  Elev <- rep(coredata$Elev[1], times=length(trees))
  averages <- aggregate(x=coredata[,10:ncol(coredata)],by=list(coredata$TreeID),FUN="mean", na.rm=T)
  nyears <- ncol(averages)-2
  colnames(averages)<- c("Tree", seq(from=yrcored, to= yrcored-nyears, by=-1))
  treedata <- cbind(Transect,Species,Elev,averages)
  
  return(treedata)
}


################################
#     read tree core data
################################
## low elevation
csv1 <- "CO-POTR-L_RWdata01_29_14.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
coredata1 <- core_clean(data=data1) #unaveraged cores (i.e. 2 cores per tree still)
treedata1 <- core_average(data=data1) # cores averaged per tree
## mid elevation
csv2 <- "CO-POTR-M_RWdata01_31_14.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
coredata2 <- core_clean(data=data2)
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "CO-POTR-H_RWdata04_14_14v2.csv"
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
coredata3 <- core_clean(data=data3)
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
rw.detH <- core_trunc(workingdata=treedata3, detmethod=method, detrend=TRUE, drop=dropyrs)
# Mid Elevation
rw.detM <- core_trunc(workingdata=treedata2, detmethod=method, detrend=TRUE, drop=dropyrs)
# Low Elevation
rw.detL <- core_trunc(workingdata=treedata1, detmethod=method, detrend=TRUE, drop=dropyrs)


### Making a .csv of master chronologies
L.master <- apply(X=rw.detL, MARGIN=1,FUN=mean, na.rm=TRUE)
M.master <- apply(X=rw.detM, MARGIN=1,FUN=mean, na.rm=TRUE)
H.master <- apply(X=rw.detH, MARGIN=1,FUN=mean, na.rm=TRUE)
POTR.master <- data.frame(L.master, M.master, H.master)
names(POTR.master) <- c("POTR.L", "POTR.M", "POTR.H")
#write.csv(POTR.master, file="MASTER_CO-POTR_04172014.csv")

### writing master chrons that I can import to WinDendro
Low<- crossdate(data=data1)
Mid <- crossdate(data=data2)
High <- crossdate(data=data3)

Ldet <- detrend(Low, method="Spline")
Lmaster <- chron(Ldet)

Mdet <- detrend(Mid, method="Spline")
Mmaster <- chron(Mdet)

Hdet <- detrend(High, method="Spline")
Hmaster <- chron(Hdet)

write.csv(Lmaster, "POTR_Lmaster_dplr.csv")
write.csv(Mmaster, "POTR_Mmaster_dplr.csv")
write.csv(Hmaster, "POTR_Hmaster_dplr.csv")
plot(Lmaster$xxxstd~as.numeric(row.names(Lmaster)), type="l", col="red")
lines(L.master~as.numeric(names(L.master)), col="darkred")
lines(Mmaster$xxxstd~as.numeric(row.names(Mmaster)), col="darkblue")

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
  #=============== data frames to use ==============
  #rw.clim.df # detrended rw and climate data by tree and year and ready to go
  rw.clim.df.narm <- na.omit(rbind(rw.clim.C, rw.clim.N)) # same df but na's removed 'casue {nlme} doesn't like them
  return(rw.clim.df.narm)
}
#=======================================================
# for dealing with scaled climate variables
Hrw.clim.df.narm <-combine_rwclim(rw.det=rw.detH,clim.des.scale=climH.des.scale,elev="H")
Mrw.clim.df.narm <-combine_rwclim(rw.det=rw.detM,clim.des.scale=climH.des.scale,elev="M")
Lrw.clim.df.narm <-combine_rwclim(rw.det=rw.detL,clim.des.scale=climH.des.scale,elev="L")


rw.clim.allelev <- rbind(Hrw.clim.df.narm, Mrw.clim.df.narm, Lrw.clim.df.narm)
#write.csv(rw.clim.allelev, file="CO-POTR-allRWIClim.csv")


# or for dealing with raw climate variables
Hrw.clim.raw <-combine_rwclim(rw.det=rw.detH,clim.des.scale=climH.des,elev="H")
Mrw.clim.raw <-combine_rwclim(rw.det=rw.detM,clim.des.scale=climH.des,elev="M")
Lrw.clim.raw <-combine_rwclim(rw.det=rw.detL,clim.des.scale=climH.des,elev="L")

rw.clim.allraw <- rbind(Hrw.clim.raw, Mrw.clim.raw, Lrw.clim.raw)
#write.csv(rw.clim.allraw, file="CO-POTR-allRWI_rawClim.csv")




#==================================================
##########     Making quick spark plots of data
#==================================================
POTR.master <- read.csv("MASTER_CO-POTR_04172014.csv")
rw.clim.allelev <- read.csv("CO-POTR-allRWIClim.csv")


#### Change this to get different elevations
toplot <- subset(x=rw.clim.allelev,subset=Elev=="L")
toplot <- toplot[order(toplot$Tree),]
# with H & L layered #THIS IS NOT WORKING, AND I DON'T KNOW WHY
sparkplots(years=toplot$Year,names=toplot$Tree, data=toplot$PPTGS, types=toplot$Comp,legend=TRUE)
# with all plotted as lines
sparkplots(years=toplot$Year,names=toplot$Tree, data=toplot$value, types=rep(x=1,times=nrow(toplot)),legend=FALSE)


#===============================================
## Making full Sparkplot figures
#===============================================
# making the plot with masters on top and parkplots of each 
# elevation below
X <- as.numeric(rownames(POTR.master))
mat <- matrix(c(1,1,1,
                2,3,4),nrow=2,ncol=3, byrow=T)

quartz(width=8, height=6)
par(mar=c(3,3,1,2), oma=c(2,3,1,1))
layout(mat=mat, heights=c(.6,1,1,1))
layout.show(4)
linecol <- "#AAAAAA"
plot(POTR.L~X, data=POTR.master, type="n", lwd=2, col="red3", ylab="", xlab="", ylim=c(-.2, 4.7), xlim=c(1900,2013), yaxt="n", mgp = c(2,.7,0), cex.axis=1.2)
#legend("bottomleft", legend=c("Low", "Mid", "High", "Individual"), lwd=c(2,2,2,1), col= c("red3", "purple3", "blue3", "#999999"), ncol=4, bty="n")
mtext(text="Detrended Master \n Chronologies",side=2, line=1.4)
text(x=c(2011,2011,2011), y=c(1.15,2.65,4.15), pos=4, labels=c("Low", "Mid", "High"), font=2)
abline(h=1, lty=2, col="gray")
abline(h=2.5, lty=2, col="gray")
abline(h=4, lty=2, col="gray")


#------ adding individual chronologies under masters -----------
toplot <- subset(x=rw.clim.allelev,subset=Elev=="L")
years <- toplot$Year
names <- toplot$Tree
data <- toplot$value
#types <- toplot$Comp
types=rep(x=1,times=nrow(toplot))
comps <- unique(toplot$Comp)

name.un <- unique(names)
name.un <- name.un[order(name.un)]
type.un <- unique(types)
year <- as.numeric(as.character(years))

for(i in 1:length(name.un)){
  n <- name.un[i] #store the name of desired tree we're looping through
  temp <- data[which(names==n)]
  Comp <- toplot$Comp[which(names==n)]
  #tempsc <- temp/max(temp)# - min(temp)) #/max(temp-min(temp)) # making values range from 0 to 1
  lines(x=year[which(names==n)],y=(temp), lwd=1, col=linecol)
}
points(POTR.L~X, data=POTR.master, type="l", lwd=2, col="red3")#, ylab="Detrended Master Chronologies", xlab="Year", ylim=c(0, 4.5), yaxt="n")
#  text (max(year), y = i-0.5, n, pos=4, xpd=NA )

toplot <- subset(x=rw.clim.allelev,subset=Elev=="M")
years <- toplot$Year
names <- toplot$Tree
data <- toplot$value
#types <- toplot$Comp
types=rep(x=1,times=nrow(toplot))
comps <- unique(toplot$Comp)

name.un <- unique(names)
name.un <- name.un[order(name.un)]
type.un <- unique(types)
year <- as.numeric(as.character(years))

for(i in 1:length(name.un)){
  n <- name.un[i] #store the name of desired tree we're looping through
  temp <- data[which(names==n)]
  Comp <- toplot$Comp[which(names==n)]
  #tempsc <- temp/max(temp)# - min(temp)) #/max(temp-min(temp)) # making values range from 0 to 1
  lines(x=year[which(names==n)],y=(temp + 1.5), lwd=1, col=linecol)
}

toplot <- subset(x=rw.clim.allelev,subset=Elev=="H")
years <- toplot$Year
names <- toplot$Tree
data <- toplot$value
#types <- toplot$Comp
types=rep(x=1,times=nrow(toplot))
comps <- unique(toplot$Comp)

name.un <- unique(names)
name.un <- name.un[order(name.un)]
type.un <- unique(types)
year <- as.numeric(as.character(years))

for(i in 1:length(name.un)){
  n <- name.un[i] #store the name of desired tree we're looping through
  temp <- data[which(names==n)]
  Comp <- toplot$Comp[which(names==n)]
  #tempsc <- temp/max(temp)# - min(temp)) #/max(temp-min(temp)) # making values range from 0 to 1
  lines(x=year[which(names==n)],y=(temp + 3), lwd=1, col=linecol)
}
points(POTR.M+1.5~X, data=POTR.master, type="l", lwd=2, col="purple3")
points(POTR.H+3~X, data=POTR.master, type="l", lwd=2, col="blue3")


# ---------- now making sparkplots below as figs 2:4 --------

par(mar=c(1,0,2,1))
toplot <- subset(x=rw.clim.allelev,subset=Elev=="L")
main <- "Low elevation"
#toplot <- toplot[order(toplot$Tree),]


years <- toplot$Year
names <- toplot$Tree
data <- toplot$value
#types <- toplot$Comp
types=rep(x=1,times=nrow(toplot))
comps <- unique(toplot$Comp)

name.un <- unique(names)
name.un <- name.un[order(name.un)]
type.un <- unique(types)
year <- as.numeric(as.character(years))


plot(1,type="n",axes=F,ann=F, xlim=c(min(year),max(year)),ylim=c(0,length(name.un)/length(type.un)), yaxs="i", xaxs="i", main=main)   #ensure axes go to the edge of the data points (remove gap between 0 and x and y axes)
axis(1)
#for(j in 1:length(type.un)){
#  t <- type.un[j]
#  namej <- names[types==t]
#  namej.un <- unique(namej)  #get the unique vector of names that are type t
for(i in 1:length(name.un)){
  n <- name.un[i] #store the name of desired tree we're looping through
  temp <- data[which(names==n)]
  Comp <- toplot$Comp[which(names==n)]
  tempsc <- temp/max(temp)# - min(temp)) #/max(temp-min(temp)) # making values range from 0 to 1
  lines(x=year[which(names==n)],y=(tempsc+i-1), lwd=1.5)
}
mtext(text=main, side=3, line=.4, cex=1, adj=.55)


#------------------ Mid elevation sparkplots

toplot <- subset(x=rw.clim.allelev,subset=Elev=="M")
main <- "Mid elevation"
years <- toplot$Year
names <- toplot$Tree
data <- toplot$value
#types <- toplot$Comp
types=rep(x=1,times=nrow(toplot))
comps <- unique(toplot$Comp)
name.un <- unique(names)
name.un <- name.un[order(name.un)]
type.un <- unique(types)
year <- as.numeric(as.character(years))
plot(1,type="n",axes=F,ann=F, xlim=c(min(year),max(year)),ylim=c(0,length(name.un)/length(type.un)), yaxs="i", xaxs="i", main=main)   #ensure axes go to the edge of the data points (remove gap between 0 and x and y axes)
axis(1)
for(i in 1:length(name.un)){
  n <- name.un[i] #store the name of desired tree we're looping through
  temp <- data[which(names==n)]
  Comp <- toplot$Comp[which(names==n)]
  tempsc <- temp/max(temp)# - min(temp)) #/max(temp-min(temp)) # making values range from 0 to 1
  lines(x=year[which(names==n)],y=(tempsc+i-1), lwd=1.5)
  #text (max(year), y = i-0.5, n, pos=4, xpd=NA )     
}
#polygon(x=c(2000.5,2000.5,2003.5,2003.5), y=c(length(name.un)/length(type.un)+.5,0,0,length(name.un)/length(type.un)+.5), col="#00000015",border="#00000015")
#polygon(x=c(1975.5,1975.5,1978.5,1978.5), y=c(length(name.un)/length(type.un)+.5,0,0,length(name.un)/length(type.un)+.5), col="#00000020",lwd=.001)
mtext(text=main, side=3, line=.4, cex=1, adj=.55)

#-------------------- High elevation sparkplots

toplot <- subset(x=rw.clim.allelev,subset=Elev=="H")
# can also knock of 11 and 12 to make each elevation have the same number of trees
#toplot <- toplot[which(toplot$Tree != "H.11H" & toplot$Tree != "H.11L" & toplot$Tree != "H.12H" & toplot$Tree != "H.12L"),]
main <- "High elevation"
years <- toplot$Year
names <- toplot$Tree
data <- toplot$value
#types <- toplot$Comp
types=rep(x=1,times=nrow(toplot))
comps <- unique(toplot$Comp)
name.un <- unique(names)
name.un <- name.un[order(name.un)]
type.un <- unique(types)
year <- as.numeric(as.character(years))
plot(1,type="n",axes=F,ann=F, xlim=c(min(year),max(year)),ylim=c(0,length(name.un)/length(type.un)), yaxs="i", xaxs="i", main=main)   #ensure axes go to the edge of the data points (remove gap between 0 and x and y axes)
axis(1)
for(i in 1:length(name.un)){
  n <- name.un[i] #store the name of desired tree we're looping through
  temp <- data[which(names==n)]
  Comp <- toplot$Comp[which(names==n)]
  tempsc <- temp/max(temp)# - min(temp)) #/max(temp-min(temp)) # making values range from min/max to 1
  lines(x=year[which(names==n)],y=(tempsc+i-1), lwd=1.5)
  #text (max(year), y = i-0.5, n, pos=4, xpd=NA )     
}
#polygon(x=c(2000.5,2000.5,2003.5,2003.5), y=c(length(name.un)/length(type.un)+.5,0,0,length(name.un)/length(type.un)+.5), col="#00000015", border="#00000015")
#polygon(x=c(1975.5,1975.5,1978.5,1978.5), y=c(length(name.un)/length(type.un)+.5,0,0,length(name.un)/length(type.un)+.5), col="#00000020",lwd=.001)
mtext(text=main, side=3, line=.4, cex=1, adj=.55)

#------------------------------------------------------------------------



#############################################################
################### Visual crossdating figure #################
#############################################################

#==============================================
# overlaying PIPO and POTR master chronologies
##==============================================
PIPO.master <- read.csv("MASTER_CO-PIPO_02072014.csv")
POTR.master <- read.csv("MASTER_CO-POTR_04172014.csv")
rw.clim.allelev <- read.csv("CO-POTR-allRWIClim.csv")


quartz(width=8, height=6)
ylim <- c(0.2, 4.7) # initially c(-0.2, 4.7)
plot(POTR.L~X, data=POTR.master, type="n", lwd=2, col="red3", ylab="", xlab="", ylim=ylim, xlim=c(1900,2013), yaxt="n", mgp = c(2,.7,0), cex.axis=1.2)
#legend("bottomleft", legend=c("Low", "Mid", "High", "Individual"), lwd=c(2,2,2,1), col= c("red3", "purple3", "blue3", "#999999"), ncol=4, bty="n")
mtext(text="Detrended Master \n Chronologies",side=2, line=1.4)
text(x=c(2011,2011,2011), y=c(1.15,2.65,4.15), pos=4, labels=c("Low", "Mid", "High"), font=2)
abline(h=1, lty=2, col="gray")
abline(h=1.5, lty=2, col="gray")
abline(h=2, lty=2, col="gray")
points(POTR.L~X, data=POTR.master, type="l", lwd=2, col="red3")
points(PIPO.L~X, data=PIPO.master, type="l", lwd=2, col="red4")
# initally +1.5
points(POTR.M+1.5~X, data=POTR.master, type="l", lwd=2, col="purple3")
points(PIPO.M+1.5~X, data=PIPO.master, type="l", lwd=2, col="purple4")
#initially +3
points(POTR.H+3~X, data=POTR.master, type="l", lwd=2, col="blue3")
points(PIPO.H+3~X, data=PIPO.master, type="l", lwd=2, col="blue4")

legend("top",legend=c("PIPO.L", "POTR.L", "PIPO.M", "POTR.M", "PIPO.H", "POTR.H"), col=c("red4","red3","purple4", "purple3", "blue4","blue3"), lwd=2, ncol=3, bty="n")

# looking at some indicator years
abline(v=2002, col="gray")
abline(v=1989, col="gray")
abline(v=1986, col="gray")
abline(v=1983, col="gray", lty=2)
abline(v=1977, col="gray")
abline(v=1969, col="gray", lty=2)
abline(v=1957, col="gray", lty=2)
abline(v=1949, col="gray")
abline(v=1947, col="gray", lty=2)
abline(v=1943, col="gray")
abline(v=1934, col="gray")
abline(v=1922, col="gray")
abline(v=1915, col="gray")
#-------------------------------------------------------



#=======================================
### Looking at individual cores
#=======================================

#plotchron function for visual crossdating
#### NOTE: I HAVE NOT automated this for different elevations. It's still hard coded
#--------------------------------------------------------------------
plotchron <- function(i, elev="H"){
  X <- as.numeric(rownames(POTR.master))
  quartz(width=8, height=6)
  par(mar=c(3,3,1,2), oma=c(2,3,1,1))
  
  #layout(mat=mat, heights=c(.6,1,1,1))
  #layout.show(4)
  linecol <- "#666666"
############ Change for different elevations
  plot(POTR.H~X, data=POTR.master, type="n", lwd=2, col="red3", ylab="", xlab="", ylim=c(-.2, 2), xlim=c(1900,2013), yaxt="n", mgp = c(2,.7,0), cex.axis=1.2, main=name.un[i])
  #legend("bottomleft", legend=c("Low", "Mid", "High", "Individual"), lwd=c(2,2,2,1), col= c("red3", "purple3", "blue3", "#999999"), ncol=4, bty="n")
  mtext(text="Detrended Master \n Chronologies",side=2, line=1.4)
  #text(x=c(2011,2011,2011), y=c(1.15,2.65,4.15), pos=4, labels=c("Low", "Mid", "High"), font=2)
  abline(h=1, lty=2, col="gray")
  #abline(h=2.5, lty=2, col="gray")
  #abline(h=4, lty=2, col="gray")
  # looking at some indicator years
  abline(v=2002, col="gray")
  abline(v=1989, col="gray")
  abline(v=1986, col="gray")
  abline(v=1983, col="gray", lty=2)
  abline(v=1977, col="gray")
  abline(v=1969, col="gray", lty=2)
  abline(v=1957, col="gray", lty=2)
  abline(v=1949, col="gray")
  abline(v=1947, col="gray", lty=2)
  abline(v=1943, col="gray")
  abline(v=1934, col="gray")
  abline(v=1922, col="gray")
  abline(v=1915, col="gray")
  
  #------ adding individual chronologies under masters -----------
  toplot <- subset(x=rw.clim.allelev,subset=Elev==elev)
  years <- toplot$Year
  names <- toplot$Tree
  data <- toplot$value
  #types <- toplot$Comp
  types=rep(x=1,times=nrow(toplot))
  comps <- unique(toplot$Comp)
  
  name.un <- unique(names)
  name.un <- name.un[order(name.un)]
  type.un <- unique(types)
  year <- as.numeric(as.character(years))
################# Change for different elevation bands  
  points(POTR.H~X, data=POTR.master, type="l", lwd=2, col="red3")
  
  
  #for(i in 1:length(name.un)){
  n <- name.un[i] #store the name of desired tree we're looping through
  temp <- data[which(names==n)]
  Comp <- toplot$Comp[which(names==n)]
  #tempsc <- temp/max(temp)# - min(temp)) #/max(temp-min(temp)) # making values range from 0 to 1
  lines(x=year[which(names==n)],y=(temp), lwd=1, col=linecol)
  #}
  points(POTR.L~X, data=POTR.master, type="l", lwd=2, col=rgb(red=.9, green=0, blue=0, alpha=.2))
}
#====================================================================


#=========================================
## actual code to plot chronologies
#=========================================
plotchron(2)

# and then throw on top the second core
i <- 1
n <- name.un[i] #store the name of desired tree we're looping through
temp <- data[which(names==n)]
Comp <- toplot$Comp[which(names==n)]
#tempsc <- temp/max(temp)# - min(temp)) #/max(temp-min(temp)) # making values range from 0 to 1
lines(x=year[which(names==n)],y=(temp), lwd=1, col="grey")


toplot <- subset(x=rw.clim.allelev,subset=Elev=="L")
# can also knock of 11 and 12 to make each elevation have the same number of trees
#toplot <- toplot[which(toplot$Tree != "H.11H" & toplot$Tree != "H.11L" & toplot$Tree != "H.12H" & toplot$Tree != "H.12L"),]
main <- "High elevation"
years <- toplot$Year
names <- toplot$Tree
data <- toplot$value
#types <- toplot$Comp
types=rep(x=1,times=nrow(toplot))
comps <- unique(toplot$Comp)
name.un <- unique(names)
name.un <- name.un[order(name.un)]
type.un <- unique(types)
year <- as.numeric(as.character(years))
plot(1,type="n",axes=F,ann=F, xlim=c(min(year),max(year)),ylim=c(0,length(name.un)/length(type.un)), yaxs="i", xaxs="i", main=main)   #ensure axes go to the edge of the data points (remove gap between 0 and x and y axes)
axis(1)
for(i in 1:length(name.un)){
  n <- name.un[i] #store the name of desired tree we're looping through
  temp <- data[which(names==n)]
  Comp <- toplot$Comp[which(names==n)]
  tempsc <- temp/max(temp)# - min(temp)) #/max(temp-min(temp)) # making values range from min/max to 1
  lines(x=year[which(names==n)],y=(tempsc+i-1), lwd=1.5)
  #text (max(year), y = i-0.5, n, pos=4, xpd=NA )     
}
#polygon(x=c(2000.5,2000.5,2003.5,2003.5), y=c(length(name.un)/length(type.un)+.5,0,0,length(name.un)/length(type.un)+.5), col="#00000015", border="#00000015")
#polygon(x=c(1975.5,1975.5,1978.5,1978.5), y=c(length(name.un)/length(type.un)+.5,0,0,length(name.un)/length(type.un)+.5), col="#00000020",lwd=.001)
mtext(text=main, side=3, line=.4, cex=1, adj=.55)











