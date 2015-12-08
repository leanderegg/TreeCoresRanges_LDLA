###########################################################
###### Plotting Climate Normals from the Three Transects
##      Normals pulled from ClimateWNA 8/15
##      Goal:produce summary figure for MS
#############################################################

# This code originally plotted wrong normals data
# that was fucked up by typos in the lat long fields
# Once I fixed those and redownloaded the normals, things changed
# also originally plotted chiricahua climate for DDIG proposal (figure never included)
# so a related version is still in that folder called FIG_Chiricahua_Climate.R or some such
# This code has been cleaned as of 11/1/15 to show updated climate Normals (and DDIG code has been deleted)


## as of 12/4/15, I've started revamping the pca stuff, but still way too many variables to interpret
# also hoping to look at 1901-1930 vs 1980-2010 normals
# -- I think I'll just have to stick with the informed winnowing of variables

require(RColorBrewer)
require(scatterplot3d)
require(rgl)
require(Rcmdr)
require(reshape2)
require(dplyr)

# optimized pairplot function
source("/Users/leeanderegg/Desktop/ZuurMixedModelling/AllRCode/HighstatLibV6.R")





###### Plotting Climate Normals from my transects ############

# relavent climate axes from Albright et al. 2013
# Annual Water Deficit vs Annual AET
# Summer Water Deficit vs April SWE
# Relevant variabes from Littel et al 2008
# July Temp vs. Jan Temp
# July Precip vs. Jan Precip
# July Temp vs. Annual Precip
# Spr-Sum Precip vs. GDD

### old versions of normal data (do not trust) ##
# climates <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WholeProjectCode/MeanLocations_Normal_1961_1990MSY.csv", header=T)
# climates2 <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WholeProjectCode/MeanLocations_Normal_1961_1990Y.csv", header=T)

### most up to date climate normals as of 11/1/15
# normals 1970-2000
climates <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WholeProjectCode/climateWNA data/Normals/forclimateWNA_MeanLocations_8_5_15_Normal_1971_2000MSY.csv", header=T)
# normals 1980-2010
climatesrecent <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WholeProjectCode/climateWNA data/Normals/forclimateWNA_MeanLocations_8_5_15_Normal_1961_1990MSY.csv", header=T)
climatesrecent$PPT_an <- with(climatesrecent, PPT_wt + PPT_sp + PPT_sm + PPT_at)
climatesrecent$DD5_an <- with(climatesrecent, DD5_wt + DD5_sp + DD5_sm + DD5_at)
climatesrecent$PPT_spsm <- with(climatesrecent, PPT_sp + PPT_sm)

# normals 1901-1930
climatesold <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WholeProjectCode/climateWNA data/Normals/forclimateWNA_MeanLocations_8_5_15_Normal_1901_1930MSY.csv", header=T)
climatesold$PPT_an <- with(climatesold, PPT_wt + PPT_sp + PPT_sm + PPT_at)
climatesold$DD5_an <- with(climatesold, DD5_wt + DD5_sp + DD5_sm + DD5_at)
climatesold$PPT_spsm <- with(climatesold, PPT_sp + PPT_sm)

# colnames(climates3)[1:2]<- c("State", "SpeciesBand")
# climates3$Species <- t(data.frame(strsplit(as.character(climates3$SpeciesBand),split = "-"))[1,])
# climates3$


###### Running PCAs on monthly and seasonal data ######
# pull out monthly climate variables and 
    # get rid of the DD_0, DD_18, and DD18
    # drop columns with all 0s (CMDs in winter)
monthlyclims <- climates[,c(8:175)]
monthlyclims <- monthlyclims[,-c(grep("DD_0", colnames(monthlyclims)), grep("DD_18", colnames(monthlyclims)), grep("DD18", colnames(monthlyclims)))]
monthlyclims <- monthlyclims[,-which(colMeans(monthlyclims)==0)]
#seasonalclims <- climates[,c(grep("_sp", colnames(climates)), grep("_sm", colnames(climates)), grep("_at", colnames(climates)),grep("_wt", colnames(climates)))]
# pull out seasonal climate variables and drop shitty DD measurements
seasonalclims <- climates[,c(176:232)]
seasonalclims <- seasonalclims[,-c(grep("DD_0", colnames(seasonalclims)), grep("DD_18", colnames(seasonalclims)), grep("DD18", colnames(seasonalclims)))]
seasonalclims <- seasonalclims[,-which(colMeans(seasonalclims)==0)]

monthlypca <- prcomp(monthlyclims,scale. = TRUE)
biplot(monthlypca) # holy shit.
summary(monthlypca) # first two PCs = 82 % of variance
seasonalpca <- prcomp(seasonalclims, scale. = T)
biplot(seasonalpca) # still tough to read
summary(seasonalpca) # first two PCs = 88% of variance
  # this no longer pulls out between mountain vs. within mountain variation. bummer.

Mypairs(seasonalclims[,c(10,11,12,13,14,15,29:35)])
## ok, this seems mostly pointless. I think most things are strongly linearaly related, often because the temp and ppt variables are used to calculate the others.

#-- Seasonal --
# Tmax_sm
# Tmin_wt
# PPT_wn PPT_sp PPT_sm PPT_at
# CMD_sp CMD_sm CMD_at

#-- Annual --
# MAT
# MAP
# DD5
# NFFD
# PAS
# CMD


col <- paste(brewer.pal (n=3, name="Set1"), "AA", sep="")  
cold <- brewer.pal (n=3, name="Set1")
coll <- paste(brewer.pal (n=3, name="Set1"), "66", sep="")
cols <- rep(col, times=xtabs(~climates$State))
colsdark <- rep(cold, times=xtabs(~climates$State))
colslight <- rep(coll, times=xtabs(~climates$State))


statepch <- c(rep(15,times=9), rep(16, times=9), rep(17, times=8))
## relavent pchs: 0,1,2,15,16,17,18
climates$pchs <- climates$Species
levels(climates$pchs) <- list("0"="ABLA", "18" = "POTR", "15"="PIPO", "16"="PSME","2" = "TSHE" )
climates$pchs <- as.numeric(as.character(climates$pchs))

#### Making plots with Littel axes #########
quartz(width=7, height=5)
par(mfrow=c(2,2), mar=c(3,3,0,0), oma=c(0,0,1,1), mgp=c(2,1,0))
# July Temp vs. Jan Temp
plot(Tave07~Tave01,data=climates, col=cols, pch=pchs)
#legend("topright",legend=c("CO", "MT", "WA"),col=col, pch=c(15,16,17), xpd=F, horiz=T)
# July Precip vs. Jan Precip
plot(PPT07~PPT01,data=climates, col=cols, pch=pchs)
# July Temp vs. Annual Precip
plot(Tave07~I(PPT_wt+PPT_sp+PPT_sm+PPT_at),data=climates, col=cols, pch=pchs, xlab="PPTan")
# Spr-Sum Precip vs. GDD
plot(I(PPT_sp + PPT_sm)~I(DD5_wt+DD5_sp+DD5_sm+DD5_at),data=climates, col=cols, pch=pchs, ylab="PPT sp+sm",xlab="annual GDD")



###### same figure as above, but with 1901-2010 differences ###
quartz(width=7, height=5, title ="1901 vs 1960")
par(mfrow=c(2,2), mar=c(3,3,0,0), oma=c(0,0,1,1), mgp=c(2,1,0))
# July Temp vs. Jan Temp
plot(Tave07~Tave01,data=climatesold, col=colslight, pch=climates$pchs)
points(Tave07~Tave01,data=climatesrecent, col=colsdark, pch=climates$pchs)
arrows(x0 = climatesold$Tave01, y0=climatesold$Tave07, x1=climatesrecent$Tave01, y1=climatesrecent$Tave07, col=cols, length=.05)
#legend("topright",legend=c("CO", "MT", "WA"),col=col, pch=c(15,16,17), xpd=F, horiz=T)
# July Precip vs. Jan Precip
plot(PPT07~PPT01,data=climatesold, col=colslight, pch=climates$pchs)
points(PPT07~PPT01,data=climatesrecent, col=colsdark, pch=climates$pchs)
arrows(x0 = climatesold$PPT01, y0=climatesold$PPT07, x1=climatesrecent$PPT01, y1=climatesrecent$PPT07, col=cols, length=.05)
# July Temp vs. Annual Precip
plot(Tave07~PPT_an,data=climatesold, col=colslight, pch=climates$pchs, xlab="annual PPT")
points(Tave07~PPT_an,data=climatesrecent, col=colsdark, pch=climates$pchs)
arrows(x0 = climatesold$PPT_an, y0=climatesold$Tave07, x1=climatesrecent$PPT_an, y1=climatesrecent$Tave07, col=cols, length=.05)

# Spr-Sum Precip vs. GDD
plot(PPT_spsm~DD5_an,data=climatesold, col=colslight, pch=climates$pchs, ylab="spring and summer PPT",xlab="annual GDD")
points(PPT_spsm~DD5_an,data=climatesrecent, col=colsdark, pch=climates$pchs)
arrows(x0 = climatesold$DD5_an, y0=climatesold$PPT_spsm, x1=climatesrecent$DD5_an, y1=climatesrecent$PPT_spsm, col=cols, length=.05)






plot(Eref~MAP, data=climates, col=cols, pch=15, cex= scale(climates$Tmax_sm-climates$Tmin_wt)+2.5)#(((climates$Tmax_sm-climates$Tmin_wt)-20.4)/16.1)*3)
legend("topleft",legend=c("CO", "MT", "WA"),col=col, pch=16, horiz=T)

plot(DD5~CMD, data=climates, col=cols, pch=sp)
legend("topright",legend=c("CO", "MT", "WA"),col=col, pch=c(15,16,17), horiz=T)


plot(PAS~CMD_sm, data=climates, col=cols, pch=sp)
legend("topright",legend=c("CO", "MT", "WA"),col=col, pch=c(15,16,17), horiz=T)

plot(Tmin_wt~CMD_sm, data=climates, col=cols, pch=sp)
legend("topright",legend=c("CO", "MT", "WA"),col=col, pch=c(15,16,17), horiz=F)

plot(Tmin_wt~CMD, data=climates, col=cols, pch=sp)

plot(Tmin_wt~Tmax_sm, data=climates, col=cols, pch=sp)
legend("topright",legend=c("CO", "MT", "WA"),col=col, pch=15, horiz=F)

plot(PPT_wt~Tmax_sm, data=climates, col=cols, pch=sp)
plot(PPT_wt~DD5, data=climates,col=cols, pch=sp)
plot(PPT_wt~CMD, data=climates, col=cols, pch=sp)
points(PPT_wt~CMD_sm, data=climates, pch=sp)

plot(PAS~DD5, data=climates, pch=sp, col=cols)
plot(PAS~CMD_sm, data=climates, pch=sp, col=cols)
plot(CMD~Eref, data=climates, pch=sp, col=cols)

## making some plots for lab group presentation
quartz(width=5, height=5)
plot(Tave01~Tave07, data=climates, pch=sp, col=colsdark, xlab="July average T", ylab="January average T")
    # strong differentiation between transects on Tave01 and within transects in Tave07
plot(PPT01~PPT07, data=climates, pch=sp, col=colsdark, xlab="July Precip", ylab="January Precip")
    # same as T, but less clean

plot(CMD_sm~NFFD, data=climates, pch=sp, col=colsdark, xlab="Growing Degree Days", ylab="summer CMD (PET-P)")
legend('topleft', legend=c("CO", "MT", "WA"), col=cold, pch=c(15,16,17), bty="n")
    # actually extremely similar climate space in these two axes --> summers warm and dry in all locations


plot(I(PPT_sp + PPT_sm + PPT_at + PPT_wt) ~ Tave07, data=climates, pch=sp, col=cols)
plot(I(PPT_sp + PPT_sm)~DD5, data=climates, pch=sp, col=cols)


with(climates, scatterplot3d(I(PPT_sp + PPT_sm + PPT_at + PPT_wt), DD5, CMD_sm, pch=pchs+14, highlight.3d = F, type="h",  color=cols,))
with(climates, scatterplot3d(I(PPT_sp + PPT_sm + PPT_at + PPT_wt), DD5, CMD_sm, pch=pchs+14, highlight.3d = F, type="h",  color=cols,))
rots <- seq(95,115, by=1)
for(i in rots){
  quartz(width=3, height=3)
  s3d <- with(climates, scatterplot3d(PPT_wt, DD5_sm, CMD_sm, pch=pchs+14, highlight.3d = F, type="h",  color=cols,angle = i, main=i))  
  s3d$points3d(x = climates$PPT_wt, y=rep(1400, times=nrow(climates)), z=climates$CMD_sm,pch=climates$pchs+14, col=colslight,cex=0.8)
  
  
}
quartz(width=5, height=5.5)
with(climates, scatterplot3d(PPT_wt, DD5_sm, CMD_sm, pch=pchs+14, highlight.3d = F, type="h",  color=cols,angle = 103))

s3d <- with(climates, scatterplot3d(PPT_wt, DD5_sm, CMD_sm, pch=pchs+14, highlight.3d = F, type="h",  color=colsdark,angle = 102, xlab="Winter Precip", ylab="Summer GDD",zlab="Summer CMD", box = T, grid=T))
s3d$points3d(x = climates$PPT_wt, y=rep(1400, times=nrow(climates)), z=climates$CMD_sm,pch=climates$pchs+14, col=colslight,cex=0.8)
s3d$points3d(x = climates$PPT_wt, y=climates$DD5_sm, z=rep(50, times=26), pch=climates$pchs+14, col=colslight, cex=0.8)
s3d$points3d(x = rep(1000, times=26), y=climates$DD5_sm, z=climates$CMD_sm, pch=climates$pchs+14, col=colslight, cex=0.8)

quartz(width=4, height=4)
plot(CMD_sm~PPT_wt, data=climates, col=colsdark, pch=pchs+14)
plot(DD5_sm~PPT_wt, data=climates, col=colsdark, pch=pchs+14)
plot(CMD_sm~DD5_sm, data=climates, col=colsdark, pch=pchs+14)






with(climates, plot3d(I(PPT_sp + PPT_sm + PPT_at + PPT_wt), DD5, CMD_sm, col=cols, pch=pchs+14, box=F))
with(climates, plot3d(DD5_sm, PPT_wt, CMD_sm, col=cols, pch=pchs+14, box=F))


with(climates, scatterplot3d(I(PPT_sp + PPT_sm + PPT_at + PPT_wt), DD5, CMD_sm, col=cols))


plot(monthlypca$x[,2]~monthlypca$x[,3], pch=sp, col=cols)
## pc1 = winter and summer growiness? ( + precip, + shoulder GDD)
# = diffs between transects
## pc2 = snow amount
# = diffs within transects
## pc3 = 

plot(seasonalpca$x[,1]~climates2$PPT_wt, pch=sp, col=cols)
seasonalpca$rotation[which(abs(seasonalpca$rotation[,3])>0.1),3]
#PC1 = warm/wet winters (-) to high evap dry summers (+) (mostly winter precip)
#PC2 = warmth of summer (mostly GDD, CMD and -PASwt/sp)
plot(seasonalpca$x[,2]~seasonalpca$x[,3], pch=sp, col=cols)
plot(PPT_wt~DD5_sm, pch=pchs+14, col=cols, climates2)






##############################################
#####    Making summary climate plot #########
#####    For CO transect for DDIG     #######
#############################################

# plan: show points in MAT~MAP space and make elipses around the 
# zones of each speces, plus have points where I sampled and scale point size by
# seasonal temp differences.
require(RColorBrewer)
col <- brewer.pal (n=3, name="Set1")   
sp <- c(rep(15,times=9), rep(16, times=9), rep(16, times=8))
COclim <- climates2[1:9,]
cols <- rep(col[c(3,1,2)], each=3)
bg <- paste(cols, "66", sep="")
Tdiff <- (COclim$Tmax_sm-COclim$Tmin_wt) #seasonal temperature differential between maxT and minT
plot(Eref~MAP, data=COclim,col=cols, bg=bg, pch=21#rep(c(17,15,16), each=3)
     , ylab="Mean Annual Temp (°C)", xlab="Mean Annual Precip (mm)"
     , cex= (Tdiff-min(Tdiff))/(max(Tdiff)-min(Tdiff))*1.9 + 1)
legend("bottomleft",legend=c("Pi. ponderosa", "Po. tremuloides", "A. lasiocarpa"),col=col,pt.bg=paste(col,"66", sep=""),pch=21,bty="n", text.font=3, pt.cex=1.5)


quartz(width=4,height=3)
par(mar=c(3.3,3.3,1,1), mgp=c(2.3,1,0))
plot(Eref~MAP, data=COclim[order(COclim$Eref),], pch=16, ylim=c(450,1050), xlim=c(400,1200)
     , ylab="PET (mm)", xlab="Mean Annual Precip (mm)"
     , cex= 1.3, type="b", lwd=0.5)
legend("bottomleft",legend=c("Pi. ponderosa", "Po. tremuloides", "A. lasiocarpa"),col=col,pt.bg=paste(col,"66", sep=""),pch=21,bty="n", text.font=3, pt.cex=1.5)





#####################################
### PCA on seasonal variables #######
#####################################
clim_season <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WholeProjectCode/MeanLocations_Normal_1961_1990Seasonal.csv", header=T)
clim_annual <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WholeProjectCode/MeanLocations_Normal_1961_1990annual.csv", header=T)

pca <- princomp(x=clim_seasonal[,8:55]) #damn. Don't have enough rows yet

pca <- princomp(x=clim_annual[,8:28])
plot(pca$loadings)

require(labdsv)
pca <- prcomp(x=clim_annual[,8:28])
summary(pca)
plot(pca)
str(pca)


test <- read.table("www.try-db.org/tmpdnld/tde2015562383.txt")
test <- read.csv("/Users/leeanderegg/Dropbox/Trade-offs project/Workbook4.csv")
test2 <- read.table("/Users/leeanderegg/Dropbox/Trade-offs project/TRY_WD_allentries.txt")







climates <- read.csv("/Users/leeanderegg/Desktop/HRL lab/DDIG 2014/ChiricahuaNormals71-00.csv", header=T)
climates <- climates[order(climates$elev),] # order them by ascending elevation
Cortez <- read.csv("/Users/leeanderegg/Desktop/HRL lab/DDIG 2014/CortezClim.csv", header=T)
BayArea <- read.csv("/Users/leeanderegg/Desktop/HRL lab/DDIG 2014/BayAreaClim_Normal61-90.csv", header=T)
# BayArea <- BayArea[,match(colnames(BayArea), colnames(climates))]
test <- rbind(climates,Cortez)
# Seasonal
# Tmax_sm
# Tmin_wt
# PPT_wn PPT_sp PPT_sm PPT_at
# CMD_sp CMD_sm CMD_at
# Annual
# MAT
# MAP
# DD5
# NFFD
# PAS
# CMD
require(RColorBrewer)
col <- brewer.pal (n=3, name="Set1")  
cols <- rep(col, times=xtabs(~climates2$State))
sp <- c(rep(15,times=9), rep(16, times=9), rep(16, times=8))


plot(Tmax_sm~Tmin_wt,data=climates2, col=cols, pch=sp)
legend("topright",legend=c("CO", "MT", "WA"),col=col, pch=15, xpd=F, horiz=T)

plot(MAT~MAP, data=climates)#, col=cols, pch=sp, cex= (((climates2$Tmax_sm-climates2$Tmin_wt)-20.4)/16.1)*3)
points(MAT~MAP, data=Cortez)
legend("topleft",legend=c("CO", "MT", "WA"),col=col, pch=15, horiz=T)

plot(DD5~CMD, data=climates2, col=cols, pch=sp)
legend("topright",legend=c("CO", "MT", "WA"),col=col, pch=15, horiz=T)


plot(PAS~CMD_sm, data=climates2, col=cols, pch=sp)
legend("topright",legend=c("CO", "MT", "WA"),col=col, pch=15, horiz=T)

plot(Tmin_wt~CMD_sm, data=climates2, col=cols, pch=sp)
legend("topright",legend=c("CO", "MT", "WA"),col=col, pch=15, horiz=F)

plot(Tmin_wt~CMD, data=climates2, col=cols, pch=sp)

plot(Tmin_wt~Tmax_sm, data=climates2, col=cols, pch=sp)
legend("topright",legend=c("CO", "MT", "WA"),col=col, pch=15, horiz=F)

plot(PPT_wt~Tmax_sm, data=climates2, col=cols, pch=sp)
plot(PPT_wt~DD5, data=climates2,col=cols, pch=sp)
plot(PPT_wt~CMD, data=climates2, col=cols, pch=sp)
points(PPT_wt~CMD_sm, data=climates2, pch=sp)



##############################################
#####    Making summary climate plot #########
#####    For Chiricahua transect for DDIG     #######
#############################################

# plan: show points in MAT~MAP space and make elipses around the 
# zones of each speces, plus have points where I sampled and scale point size by
# seasonal temp differences.
require(RColorBrewer)
quartz(width=4.5, height=3.5)
par(mar=c(4.5,5,1,1))
c <- brewer.pal (n=7, name="RdBu")   # get the two red and blue end colors that I like
col <- c[c(7,1)] # get rid of the green color
#col <- c("#EE0000", "#0000EE") #just for diagnosing problems
ramp <- colorRampPalette(col) #make a function that interpolates between them
#sp <- c(rep(15,times=9), rep(16, times=9), rep(16, times=8))
#cols <- rep(col[c(3,1,2)], times=xtabs(~COclim$Species)[1:3])
#bg <- paste(cols, "66", sep="")
Tdiff <- (climates$Tmax_sm-climates$Tmin_wt) #seasonal temperature differential between maxT and minT
TdBA <- BayArea$Tmax_sm-BayArea$Tmin_wt
#Tdiff <- (test$Tmax_sm-test$Tmin_wt)
## interpolate 100 steps, round the Tdiff into 100s and then pull the propper color
colors <- ramp(101)[round((Tdiff-min(Tdiff))/(max(Tdiff-min(Tdiff)))*100)+1] 
colors <- with(climates, ramp(101)[round((Tmax_sm-min(Tmax_sm))/max(Tmax_sm-min(Tmax_sm))*100)+1])
bg <- paste(colors, "CC", sep="")
plot(MAT~MAP, data=climates,col=colors, bg=bg, pch=21#rep(c(17,15,16), each=3)
     , ylab="Mean Annual Temp (°C)", xlab="Mean Annual Precip (mm)"
     , cex= (Tdiff-min(Tdiff))/(max(Tdiff)-min(Tdiff))*2.1 + 1.2
     , ylim= c(7,16), xlim=c(400,1150))
points(x=c(480,570,640), y=c(8.6,8.4,8.2), cex=c(.98,.5,0)*2.1+1.2, pch=21, col=rev(colors[c(7,4,1)]), bg = rev(bg[c(7,4,1)]))
text(x=c(480,560,640),y=c(7.6,7.6,7.6),labels=c("35°", "32.5°", "30°"), cex=.8)
text(x=560, y = 10.2, labels="Seasonal Temp Variation", cex=.7)
text(x=560, y = 9.5, labels="(Summer Tmax - Winter Tmin)", cex=.6)
points(MAT~MAP, data=BayArea[1,], pch='S')#cex=1.2+(TdBA-min(Tdiff))/(max(TdBA)-min(TdBA))*2.1 + 1.2)
points(MAT~MAP, data=BayArea[2,], pch='B')
points(MAT~MAP, data=BayArea[3,], pch='R')
points(MAT~MAP, data=BayArea[4,], pch='E')

#plot(MAT~MAP, data=BayArea, cex=(TdBA-min(TdBA))/(max(Tdiff)-min(Tdiff))*2.1 + 1.2)



# Adding in California and CO
##############################################
#####    Making summary climate plot #########
#####    For Chiricahua transect for DDIG     #######
#############################################

# plan: show points in MAT~MAP space and make elipses around the 
# zones of each speces, plus have points where I sampled and scale point size by
# seasonal temp differences.
require(RColorBrewer)
quartz(width=4.5, height=3.5)
par(mar=c(4.5,5,1,1))
c <- brewer.pal (n=7, name="RdBu")   # get the two red and blue end colors that I like
col <- c[c(7,1)] # get rid of the green color
#col <- c("#EE0000", "#0000EE") #just for diagnosing problems
ramp <- colorRampPalette(col) #make a function that interpolates between them
#sp <- c(rep(15,times=9), rep(16, times=9), rep(16, times=8))
#cols <- rep(col[c(3,1,2)], times=xtabs(~COclim$Species)[1:3])
#bg <- paste(cols, "66", sep="")
Tdiff <- (climates$Tmax_sm-climates$Tmin_wt) #seasonal temperature differential between maxT and minT
TdBA <- BayArea$Tmax_sm-BayArea$Tmin_wt
#Tdiff <- (test$Tmax_sm-test$Tmin_wt)
## interpolate 100 steps, round the Tdiff into 100s and then pull the propper color
colors <- ramp(101)[round((Tdiff-min(Tdiff))/(max(Tdiff-min(Tdiff)))*100)+1] 
colors <- with(climates, ramp(101)[round((Tmax_sm-min(Tmax_sm))/max(Tmax_sm-min(Tmax_sm))*100)+1])
bg <- paste(colors, "CC", sep="")
plot(MAT~MAP, data=climates,col=colors, bg=bg, pch=21#rep(c(17,15,16), each=3)
     , ylab="Mean Annual Temp (°C)", xlab="Mean Annual Precip (mm)"
     , cex= (Tdiff-min(Tdiff))/(max(Tdiff)-min(Tdiff))*2.1 + 1.2
     , ylim= c(7,16), xlim=c(350,1150))
points(x=c(480,570,640), y=c(8.6,8.4,8.2), cex=c(.98,.5,0)*2.1+1.2, pch=21, col=rev(colors[c(7,4,1)]), bg = rev(bg[c(7,4,1)]))
text(x=c(480,560,640),y=c(7.6,7.6,7.6),labels=c("35°", "32.5°", "30°"), cex=.8)
text(x=560, y = 10.2, labels="Seasonal Temp Variation", cex=.7)
text(x=560, y = 9.5, labels="(Summer Tmax - Winter Tmin)", cex=.6)
points(MAT~MAP, data=BayArea[1,], pch='S')#cex=1.2+(TdBA-min(Tdiff))/(max(TdBA)-min(TdBA))*2.1 + 1.2)
points(MAT~MAP, data=BayArea[2,], pch='B')
points(MAT~MAP, data=BayArea[3,], pch='R')
points(MAT~MAP, data=BayArea[4,], pch='E')
points(MAT~MAP, data=Cortez,pch="C")
seattle <- c(11.4, 927, 20.6, 360)
points(seattle[1]~seattle[2], pch="W")
#plot(MAT~MAP, data=BayArea, cex=(TdBA-min(TdBA))/(max(Tdiff)-min(Tdiff))*2.1 + 1.2)



### same plot, TMAX and CMD as axes
plot(Tmax_sm~CMD, data=climates,col=colors, bg=bg, pch=21#rep(c(17,15,16), each=3)
     , ylab="Max Summer T (°C)", xlab="CMD (mm)"
     , cex= (Tdiff-min(Tdiff))/(max(Tdiff)-min(Tdiff))*2.1 + 1.2
     , ylim= c(22,34), xlim=c(250,1150))
points(x=c(480,570,640), y=c(8.6,8.4,8.2), cex=c(.98,.5,0)*2.1+1.2, pch=21, col=rev(colors[c(7,4,1)]), bg = rev(bg[c(7,4,1)]))
text(x=c(480,560,640),y=c(7.6,7.6,7.6),labels=c("35°", "32.5°", "30°"), cex=.8)
text(x=560, y = 10.2, labels="Seasonal Temp Variation", cex=.7)
text(x=560, y = 9.5, labels="(Summer Tmax - Winter Tmin)", cex=.6)
points(Tmax_sm~CMD, data=BayArea[1,], pch='S')#cex=1.2+(TdBA-min(Tdiff))/(max(TdBA)-min(TdBA))*2.1 + 1.2)
points(Tmax_sm~CMD, data=BayArea[2,], pch='B')
points(Tmax_sm~CMD, data=BayArea[3,], pch='R')
points(Tmax_sm~CMD, data=BayArea[4,], pch='E')
points(Tmax_sm~CMD, data=Cortez,pch="C")
seattle <- c(11.4, 927, 20.6, 23, 360)
points(seattle[4]~seattle[5], pch="W")
