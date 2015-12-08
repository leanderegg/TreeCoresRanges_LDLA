###################################################
#       Plotting Relative Abundance
#           and then Regen 
#                For all 3 transects
###################################################
# LDL Anderegg
# 4/25/14
# code from CO_assumptions check, plus WA-AbundanceFigure, plus MT-AbundanceFigure

### What I did: cleaned up code from the various figures in the above files
# then set the final abundance objects to different things so I could keep them all at once
# then plotted the three transects in a three panel figure.

# also now includes a figure with regen panels for each transect

require(vegan)
require(RColorBrewer)
require(reshape)

#===================================================
#           Creating Data frames to work from
#===================================================

#                   COLORADO
#====================================================
abun <- read.csv(file="/Users/leeanderegg/Dropbox/treecoresranges(leander) (1)/CO_RingWidths/CO_Checks/CO_abundancetransects3_7.csv", header=T)

#limiting down to only live trees
abunlive <- abun[,-grep("d", x=names(abun), ignore.case=FALSE)]


# then turn into relative abundance
abunrel <- decostand(x=as.matrix(abunlive[,-1]), method="total", margin=1)
elevfloor <- round(x=abunlive$elevation*2, digits=-2)/2

# get mean relative abundance of each species for each elevation band
mabunrel <- by(abunrel,elevfloor, colMeans)

mAbunRel <- mabunrel[[1]]
for (i in 2:dim(mabunrel)){
  mAbunRel <- cbind(mAbunRel, mabunrel[[i]]) 
}
colnames(mAbunRel) <- names(mabunrel)


### melting full abundance dataset for loess smoothing
tmp <- t(mAbunRel)  
mAbunRel.df <- melt (tmp,variable_name=colnames(tmp))
names(mAbunRel.df) <- c("Elev", "Species", "RelAbun")

### Save the CO data
CmAbunRel.df2 <- mAbunRel.df[which(mAbunRel.df$Species!="PSME"),]
levels(CmAbunRel.df2$Species) <- list(PIPO = "PIPO", POTR="POTR", ABLA="ABLA", PIEN = "Picea_sp", PIED = "PIED", JUNI = "Juniperus_sp")



#                  MONTANTA
#===================================================
abun <- read.csv(file="/Users/leeanderegg/Dropbox/treecoresranges(leander) (1)/MT_RingWidths/MT_Checks/MT_abundancetransects3-28.csv", header=T)

#limiting down to only live trees
abunlive <- abun[,-grep("d", x=names(abun), ignore.case=FALSE)]
abunlive <- abunlive[,-c(8,10)] # get rid of ABGR and PIEN (and last notes column)
require(vegan) #to get decostand() standardization function

# then turn into relative abundance
abunrel <- decostand(x=as.matrix(abunlive[,-1]), method="total", margin=1)
elevfloor <- round(x=abunlive$elevation*2, digits=-2)/2

# get mean relative abundance of each species for each elevation band
mabunrel <- by(abunrel,elevfloor, colMeans)

mAbunRel <- mabunrel[[1]]
for (i in 2:dim(mabunrel)){
  mAbunRel <- cbind(mAbunRel, mabunrel[[i]]) 
}
colnames(mAbunRel) <- names(mabunrel)

### melting full abundance dataset for loess smoothing
tmp <- t(mAbunRel)  
require(reshape)
MmAbunRel.df <- melt (tmp,variable_name=colnames(tmp))
names(MmAbunRel.df) <- c("Elev", "Species", "RelAbun")




#                WASHINGTON
#===================================================
abun <- read.csv(file="/Users/leeanderegg/Dropbox/treecoresranges(leander) (1)/WA_RingWidths/WA_Checks/WA_abundancetransects3-30.csv", header=T)

#limiting down to only live trees
abunlive <- abun[,-grep("d", x=names(abun), ignore.case=FALSE)]

require(vegan) #to get decostand() standardization function

# then turn into relative abundance
abunrel <- decostand(x=as.matrix(abunlive[,-1]), method="total", margin=1)
#elevfloor <- round(x=abunlive$elevation*2, digits=-2)/2
# using my hand rounded elevations from the spreadsheet
elevfloor <- abun$Elev
# get mean relative abundance of each species for each elevation band
# note, this is getting mean relative abundance of all transects for each species, not adding all together then relativizing
mabunrel <- by(abunrel,elevfloor, colMeans)

mAbunRel <- mabunrel[[1]]
for (i in 2:dim(mabunrel)){
  mAbunRel <- cbind(mAbunRel, mabunrel[[i]]) 
}
colnames(mAbunRel) <- names(mabunrel)


### melting full abundance dataset for loess smoothing
tmp <- t(mAbunRel)  
require(reshape)
WmAbunRel.df <- melt (tmp,variable_name=colnames(tmp))
names(WmAbunRel.df) <- c("Elev", "Species", "RelAbun")




#===================================================
#           Plotting three transects
#===================================================


quartz(width=6.5, height=6.5)
par(mar=c(3,5,2,0), oma=c(2,0,2,2), mfrow=c(3,1))

cex.axis <- 1.1
cex.names <-1.1
cex.lab <- .8


#           Plotting CO transect
#===================================================
cols <- brewer.pal(n=6, name="Set1")
cols[1] <- brewer.pal(n=9, name="PuRd")[7] #changing pondi color
cols[2] <- brewer.pal(n=9, name="PuBu")[9]
cols[6] <- brewer.pal(n=8, name="Accent")[8] # change yellow to grey
span <- 0.5
plot(RelAbun~Elev, data=CmAbunRel.df2, type="n", ylab="Relative Abundance", yaxs="i", ylim=c(0,1.1), xlim=c(2200, 3550), bty="l", cex.axis=cex.axis, main="Colorado")
# can't actually get the species in the right order because of unique. need to change this up
#mtext(side=2, text="Relative Abundance", line=3.1)
species <- unique(CmAbunRel.df2$Species)[c(4,3,1,2,5,6)]
for (i in 1:length(species))
{
  sp <- species[i]
  test <- loess(RelAbun~Elev, data=CmAbunRel.df2[CmAbunRel.df2$Species==sp,],span=span)
  lines(test$fitted ~ CmAbunRel.df2$Elev[CmAbunRel.df2$Species==sp], lwd=3, col=cols[i])
}
#par(font=3)
#legend("topright",legend=c("P. ponderosa", "P. tremuloides", "A. lasiocarpa", "Picea spp.", "P. menziesii", "Juniperus sp.", "P. edulis"),col=cols, lwd=3, bty="n", cex=.9,ncol=2)
#par(font=1)
text(x=c(2410,2810, 3300,3030, 2450, 2380), y=c(0.7,0.7, 0.3, 0.3, 0.12,0.28), labels=c("*Pinus\nponderosa","*Populus\ntremuloides", "*Abies\nlasiocarpa", "Picea spp.", "P. edulis", "J. osteosperma"),col=cols, font=c(4,4,4,3,3,3), cex=c(rep(cex.names, 3), rep(cex.lab,3)))
text(x=c(2200,3550), y=c(0.5,0.4), labels = c("Pinion/Ponderosa\ntransition", "Tree Line"), srt=90 )




#           Plotting MT transect
#===================================================
span <- 0.6
cols <- brewer.pal(n=7, name="Set1")
cols[6] <- brewer.pal(n=8, name="Accent")[8]
cols[7] <- brewer.pal(n=8, name="Set1")[8]
cols <- cols[c(1,4,5,6,2,3,7)]
plot(RelAbun~Elev, data=MmAbunRel.df, type="n", ylab="Relative Abundance", yaxs="i", ylim=c(0,1.1), xlim=c(950, 2250), mgp = c(3,1,0), bty="l", cex.axis=cex.axis, main="Montana")
for (i in 1:length(unique(MmAbunRel.df$Species)))
{
  sp <- unique(MmAbunRel.df$Species)[i]
  test <- loess(RelAbun~Elev, data=MmAbunRel.df[MmAbunRel.df$Species==sp,],span=span, degree=2)
  lines(test$fitted ~ MmAbunRel.df$Elev[MmAbunRel.df$Species==sp], lwd=3, col=cols[i])
}
par(font=3)
#legend("topright", legend =unique(mAbunRel.df$Species), col=cols, lwd=3, bty="n" )
legend(x=1850, y=0.5,legend=c("Thuja plicata", "Larix occidentalis", "Pinus spp.", "Picea engelmannii"),col=cols[c(2,3,4,7)], lwd=3, bty="n", cex=cex.lab)
par(font=1)
#locs <- locator(n=3)
locsx <- c(1030, 1317, 1583)
locsy <- c(0.604, 0.838, 0.893)
#locs$x[2] <- 1270
text(x=locsx, y=locsy, labels=c("*Tsuga\nheterophylla", "*Pseudotsuga\nmenziesii", "*Abies\nlasiocarpa"),col=cols[c(1,5,6)], font=4,pos=4, cex=cex.names)
text(x=c(950,2250), y=c(0.4,0.4), labels = c("Lake McDonald", "Tree Line"), srt=90 )



#===================================================
#           Plotting WA Abundance

span <- 0.7
cols <- brewer.pal(n=7, name="Set1")
cols[6] <- brewer.pal(n=8, name="Accent")[8]
cols[c(2,3,4)] <- cols[c(4,2,3)]
plot(RelAbun~Elev, data=WmAbunRel.df, type="n", ylab="Mean Relative Abundance", xlab="", yaxs="i", ylim=c(0,1.1), xlim=c(150,1550), mgp = c(2.5,1,0), bty="l", cex.axis=cex.axis, main="Washington")
for (i in 1:length(unique(WmAbunRel.df$Species)))
{
  sp <- unique(WmAbunRel.df$Species)[i]
  test <- loess(RelAbun~Elev, data=WmAbunRel.df[WmAbunRel.df$Species==sp,],span=span)
  lines(test$fitted ~ WmAbunRel.df$Elev[WmAbunRel.df$Species==sp], lwd=3, col=cols[i])
}
par(font=3)
#legend("topleft", legend =unique(mAbunRel.df$Species), col=cols, lwd=3, bty="n" )
#legend("topright",legend=c("P. ponderosa", "P. termuloides", "A. lasiocarpa", "Picea sp.", "P. edulis", "Juniperus sp.", "P. menziesii"),col=cols[c(5,3,1,2,7,6,4)], lwd=3, bty="n", cex=.9)
par(font=1)

xlocs <- c(210.6710,  384.9774,  629.4707, 1164.3944,  863.4570, 1379.1790, 330.6715)
ylocs <- c(0.5222512, 0.2242332, 1.0154304, 0.7749657, 0.3172240, 0.1679085, 0.098828)
text(x=xlocs, y=ylocs, labels=c("*Tsuga\nheterophylla","Thuja plicata", "*Pseudotsuga menziesii", "*Abies\nlasiocarpa", "Abies grandis", "Abies\namabilis", "Acer\n macrophyllum"),col=cols, font=c(4,3,4,4,3,3,3),pos=4, cex=c(cex.names,cex.lab,cex.names,cex.names,cex.lab,cex.lab,cex.lab))
text(x=c(150,1550), y=c(0.4,0.4), labels = c("Elwha River", "Tree Line"), srt=90 )

#==============================================================

mtext(side=1, text="Elevation (m)", line=3)







#========================================================
###################### Regen Plot #######################
#========================================================

# CO data (from CO_assumptionchecks.R)
COtreesnew <- read.csv("/Users/leeanderegg/Dropbox/treecoresranges(leander) (1)/CO_RingWidths/CO_Checks/CO_focaltree_3_6.csv", header=T)
levels(COtreesnew$Band) <- list (H = "H", L = "L", M= "M") # had an empty factor in there somehow that was fucking me up
# designate appropriate species
spL <- "PIPO"
spM <- "POTR"
spH <- "ABLA"
#first going to break things down to species for easy handling
spLall <- subset(COtreesnew,subset=Species==spL)
#which(names(spLall)=="PIPO_r1")
spMall <- subset(COtreesnew,subset=Species==spM)
#which(names(spM)=="POTR_r1")
spHall <- subset(COtreesnew,subset=Species==spH)
grep(pattern=paste(spH,"_r", sep=""),x=names(spHall))

#making a df with just regen, dropping dead
spLreg <- spLall[,c(2:6,grep(pattern=paste(spL,"_r", sep=""),x=names(spLall)))]
spMreg <- spMall[,c(2:6,grep(pattern=paste(spM,"_r", sep=""),x=names(spMall)))]
spHreg <- spHall[,c(2:6,grep(pattern=paste(spH,"_r", sep=""),x=names(spHall)))]
## Indexing note: columns 6:9 are the living regen columns
columns <- 6:9 # insert index for desired columns here

regLtmp <- by(spLreg[,columns], spLreg$Band, colMeans) # get the mean per plot regen
mspLreg <- rbind(regLtmp[[2]], regLtmp[[3]], regLtmp[[1]]) # put it into rows = L, M, Hcolumns=regen classes
regMtmp <- by(spMreg[,6:9], spMreg$Band, colMeans)
mspMreg <- rbind(regMtmp[[2]], regMtmp[[3]], regMtmp[[1]])
regHtmp <- by(spHreg[,6:9], spHreg$Band, colMeans, na.rm=TRUE)
mspHreg <- rbind(regHtmp[[2]], regHtmp[[3]], regHtmp[[1]])

COmspLreg <- mspLreg / 5^2*pi * 10000  # changing to saplings/ ha
COmspMreg <- mspMreg / 5^2*pi * 10000  # in saplings/plot / m2/plot * m2/ha
COmspHreg <- mspHreg / 5^2*pi * 10000 



# MT regen (from MT-Checks.R)
MTtreesnew <- read.csv("/Users/leeanderegg/Dropbox/treecoresranges(leander) (1)/MT_RingWidths/MT_Checks/MT_focaltree_4-5.csv", header=T)
# designate species
spL <- "TSHE"
spM <- "PSME"
spH <- "ABLA"
#first going to break things down to species for easy handling
spLall <- subset(MTtreesnew,subset=Species==spL)
#which(names(spLall)=="PIPO_r1")
spMall <- subset(MTtreesnew,subset=Species==spM)
#which(names(spM)=="POTR_r1")
spHall <- subset(MTtreesnew,subset=Species==spH)
grep(pattern=paste(spH,"_r", sep=""),x=names(spHall))

#making a df with just regen, dropping dead
spLreg <- spLall[,c(2:6,grep(pattern=paste(spL,"_r", sep=""),x=names(spLall)))]
spMreg <- spMall[,c(2:6,grep(pattern=paste(spM,"_r", sep=""),x=names(spMall)))]
# going to drop out size class t
spMreg <- spMreg[,-10]
spHreg <- spHall[,c(2:6,grep(pattern=paste(spH,"_r", sep=""),x=names(spHall)))]
## Indexing note: columns 6:10 are the living regen columns
columns <- 6:9 # insert index for desired columns here
# going to run this without size class 5 first

regLtmp <- by(spLreg[,columns], spLreg$Band, colMeans) # get the mean per plot regen
mspLreg <- rbind(regLtmp[[2]], regLtmp[[3]], regLtmp[[1]]) # put it into rows = L, M, Hcolumns=regen classes
regMtmp <- by(spMreg[,6:9], spMreg$Band, colMeans)
mspMreg <- rbind(regMtmp[[2]], regMtmp[[3]], regMtmp[[1]])
regHtmp <- by(spHreg[,6:9], spHreg$Band, colMeans, na.rm=TRUE)
mspHreg <- rbind(regHtmp[[2]], regHtmp[[3]], regHtmp[[1]])

MTmspLreg <- mspLreg / 5^2*pi * 10000  # changing to saplings/ ha
MTmspMreg <- mspMreg / 5^2*pi * 10000  # in saplings/plot / m2/plot * m2/ha
MTmspHreg <- mspHreg / 5^2*pi * 10000 


# WA data (from WA-checks.R)
WAtreesnew <- read.csv("/Users/leeanderegg/Dropbox/treecoresranges(leander) (1)/WA_RingWidths/WA_Checks/WA_focaltree_4-18.csv", header=T)
levels(MTtreesnew$Band) <- list (H = "H", L = "L", M= "M") # had an empty factor in there somehow that was fucking me up
#designate species
spL <- "TSHE"
spM <- "PSME"
spH <- "ABLA"
#first going to break things down to species for easy handling
spLall <- subset(WAtreesnew,subset=Species==spL)
#which(names(spLall)=="PIPO_r1")
spMall <- subset(WAtreesnew,subset=Species==spM)
#which(names(spM)=="POTR_r1")
spHall <- subset(WAtreesnew,subset=Species==spH)
names(spHall)[grep(pattern=paste(spH,"_r", sep=""),x=names(spHall))]

#making a df with just regen, dropping dead
spLreg <- spLall[,c(2:6,grep(pattern=paste(spL,"_r", sep=""),x=names(spLall)))]
spMreg <- spMall[,c(2:6,grep(pattern=paste(spM,"_r", sep=""),x=names(spMall)))]
# going to drop out size class t
spMreg <- spMreg[,-10]
spHreg <- spHall[,c(2:6,grep(pattern=paste(spH,"_r", sep=""),x=names(spHall)))]
## Indexing note: columns 6:10 are the living regen columns
columns <- 6:9 # insert index for desired columns here
# going to run this without size class 5 first

# this was originally in the wrong order because somewhere in the script I pulled from I must have changed the levles
regLtmp <- by(spLreg[,columns], spLreg$Band, colMeans) # get the mean per plot regen
mspLreg <- rbind(rep(NA, time=4), regLtmp[[3]], regLtmp[[1]]) # put it into rows = L, M, Hcolumns=regen classes
regMtmp <- by(spMreg[,6:9], spMreg$Band, colMeans)
mspMreg <- rbind(regMtmp[[2]], regMtmp[[3]], regMtmp[[1]])
regHtmp <- by(spHreg[,6:9], spHreg$Band, colMeans, na.rm=TRUE)
mspHreg <- rbind(regHtmp[[2]], regHtmp[[3]], regHtmp[[1]])

WAmspLreg <- mspLreg / 5^2*pi * 10000  # changing to saplings/ ha
WAmspMreg <- mspMreg / 5^2*pi * 10000  # in saplings/plot / m2/plot * m2/ha
WAmspHreg <- mspHreg / 5^2*pi * 10000 



#_____________ Plotting the Figure _____________________________



quartz(width=6.5, height=6.7)
mat <- matrix(c(1,2,3,
                4,5,6,
                7,8,9),nrow=3,ncol=3, byrow=T) 
par(mar=c(3,2,1,0), oma=c(2,4,3,2))

layout(mat=mat) #,heights=c(.8,.8,.8, .8, .8,.8,1))
#layout.show(9)  #shows what the figures borders are 


cex.axis <- 1.1
cex.names <-1.1
cex.lab <- .8

#================= Plots 1:3 = CO =============
cols <- brewer.pal(n=6, name="Set1")
cols[1] <- brewer.pal(n=9, name="PuRd")[7] #changing pondi color
cols[2] <- brewer.pal(n=9, name="PuBu")[9]
cols[6] <- brewer.pal(n=8, name="Accent")[8] # change yellow to grey

names <- c("Low", "Mid", "High")
ylims <- c(0, 31000) # 26000
ylabel <- "Saplings per ha"
colstmp = brewer.pal(4,"RdPu")
barplot(t(COmspLreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[1])#, legend.text = c("size 1", "size 2", "size 3", "size 4"))
mtext(side=2, text=ylabel,line=2.5, cex=.8)
mtext(text="Colorado", side=2, line=4, font=2, cex=1.2)
text(x=7, y=28000,labels="Pi. ponderosa", font=3,pos=1, cex=1.2)
#legend("topright", legend=c("<0.5m height", "<1.3m height", "<2.5cm DBH", "2.5-5cm DBH"), fill=colstmp, border=cols[1], bty="n", cex=.8)
colstmp = brewer.pal(4, "Purples")
barplot(t(COmspMreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[2], yaxt="n")
text(x=7, y=29000,labels="Po. tremuloides", font=3,pos=1, cex=1.2)
#text(x= 8, y=37000, labels="Colorado", xpd=NA, font=2, cex=1.5, pos=1)
colstmp= brewer.pal(4, "Greens")
barplot(t(COmspHreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[3], yaxt="n")
text(x=7, y=28000,labels="A. lasiocarpa", font=3,pos=1, cex=1.2)
#text(x= 17, y=34000, labels="Colorado", xpd=NA, font=2, cex=1.5, pos=2)


#================= Plots 4:6 = MT =============
cols <- brewer.pal(n=6, name="Set1")
cols[6] <- brewer.pal(n=8, name="Accent")[8]

names <- c("Low", "Mid", "High")
ylims <- c(0, 31000)
ylabel <- "Saplings per ha"
colstmp = brewer.pal(4,"Reds")
barplot(t(MTmspLreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[1])#, legend.text = c("size 1", "size 2", "size 3", "size 4"))
mtext(side=2, text=ylabel,line=2.5, cex=.8)
mtext(text="Montana", side=2, line=4, font=2, cex=1.25)
text(x=7, y=28000,labels="T. heterophylla", font=3,pos=1, cex=1.2)
#legend("topright", legend=c("<0.5m height", "<1.3m height", "<2.5cm DBH", "2.5-5cm DBH"), fill=colstmp, border=cols[1], bty="n", cex=.8)
colstmp = brewer.pal(4, "Blues")
barplot(t(MTmspMreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[2], yaxt="n")
text(x=7, y=28000,labels="P. menziesii", font=3,pos=1, cex=1.2)
#text(x= 8, y=35000, labels="Montana", xpd=NA, font=2, cex=1.5)
colstmp= brewer.pal(4, "Greens")
barplot(t(MTmspHreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[3], yaxt="n")
#text(x= 17, y=35000, labels="Montana", xpd=NA, font=2, cex=1.5, pos=2)


#================= Plots 7:9 = WA =============
cols <- brewer.pal(n=6, name="Set1")
cols[6] <- brewer.pal(n=8, name="Accent")[8]

names <- c("Low", "Mid", "High")
ylims <- c(0, 31000)
ylabel <- "Saplings per ha"
colstmp = brewer.pal(4,"Reds")
barplot(t(WAmspLreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[1])#, legend.text = c("size 1", "size 2", "size 3", "size 4"))
mtext(side=2, text=ylabel,line=2.5, cex=.8)
mtext(text="Washington", side=2, line=4, font=2, cex=1.25)
legend("left", legend=c("<0.5m height", "<1.3m height", "<2.5cm DBH", "2.5-5cm DBH"), fill=colstmp, border=cols[1], bty="n", cex=.8)
colstmp = brewer.pal(4, "Blues")
barplot(t(WAmspMreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[2], yaxt="n")
#text(x= 7, y=35000, labels="Washington", xpd=NA, font=2, cex=1.5)
colstmp= brewer.pal(4, "Greens")
barplot(t(WAmspHreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[3], yaxt="n")
#text(x= 17, y=35000, labels="Washington", xpd=NA, font=2, cex=1.5, pos=2)

