#####################################
####### Testing WA assumptions ######
############# 4/18/14 ##############
#(but mostly making mort and regen plots)


WAtreesnew <- read.csv("/Users/leeanderegg/Dropbox/treecoresranges(leander) (1)/WA_RingWidths/WA_Checks/WA_focaltree_4-18.csv", header=T)
# xtabs(~Species + Band + Comp, WAtreesnew)
# good. all trees seem to be there.



#=========== Using new, clean, and complete Dataz ===================
boxplot(DBH~Band+Species, data=WAtreesnew)
# supports what I was thinking for the DBH trends above
boxplot(BA_tot~Comp+Band+Species, data=WAtreesnew, col=c("gray", "white"))

levels(WAtreesnew$Band) <- list(L = "L", M="M", H="H") # reordering the levels so they plot low to high
levels(WAtreesnew$Species) <- list(TSHE="TSHE", PSME="PSME", ABLA="ABLA")
quartz(width=7, height=5)
par(mfrow=c(1,3),mar=c(5,5,3,0))
for(i in 1:3){
  sp <- levels(WAtreesnew$Species)[i]
  boxplot(BA_tot~Comp + Band, data=WAtreesnew[which(WAtreesnew$Species==sp),], main= sp, col=c("grey", "white"), ylim=c(0,600), ylab=expression(paste("Stand Basal Area (", ft^2, ha^-1, ")")))
}

## code for getting a histogram of BA_tot comp differences for a single band
tmp <- WAtreesnew[which(WAtreesnew$Species=="PSME" & WAtreesnew$Band=="H"),]
difs <- tmp$BA_tot[seq(from=1,to=nrow(tmp),by=2)] - tmp$BA_tot[seq(from=2,to=nrow(tmp),by=2)]
hist (difs)
  # it looks like I somehow totally shanked it for PSME-H.

bandtest <- WAtreesnew$Band
levels(bandtest)<- list(L = "L", M="M", H="H")
# ok, so high and low competition comparisons are actually pretty straight forward here
  # except for high elevation PSME. 

#===============================================
#====== checking out mortality roughly =========
quartz(width=7, height=5)
par(mfrow=c(1,3),mar=c(5,5,3,0))
for(i in 1:3){
  sp <- levels(WAtreesnew$Species)[i]
  boxplot(BA_totR~ Band, data=WAtreesnew[which(WAtreesnew$Species==sp),], 
          main= sp, col=c("grey", "white"), 
          ylab=expression(paste("Stand Basal Area (", ft^2, ha^-1, ")")))
}

## Making matrix of average % BA mortality for each species & Band
mort <- WAtreesnew$BA_sameRG/WAtreesnew$BA_same_all # getting percentage dead of total same basal area
mort[which(is.na(mort))] <- 0
perc_mort <- tapply(mort, INDEX=list(WAtreesnew$Species, WAtreesnew$Band), FUN=mean)

#plotting
quartz(width=7, height=5)
par(mfrow=c(1,3), oma=c(0,2,1,2), mar=c(5,3,3,0))
barplot (perc_mort[1,], ylim=c(0,0.30), space=1)
mtext(side=3, text=rownames(perc_mort)[1])
mtext(side=2, text= "% basal area dead", padj=-2.7)
barplot(perc_mort[2,], ylim=c(0, 0.3), yaxt="n", space=1)
mtext(side=3, text=rownames(perc_mort)[2])
barplot(perc_mort[3,], ylim=c(0, 0.3), yaxt="n", space=1)
mtext(side=3, text=rownames(perc_mort)[3])


#### Making barplot of mean BA with mortality colored in ###
BAsame <- tapply(WAtreesnew$BA_same_all, INDEX=list(WAtreesnew$Species, WAtreesnew$Band), FUN=mean)
#BAsame <- BAsame[,c(2,3,1)] # reorder elevations to go Low Med High
#BAsame <- BAsame[c(2,3,1),]
BAln <- log(BAsame)

BAdead <- BAsame * perc_mort # and make a vector that's the average dead BA

#-=--------- Focal species BA and BA dead ------------
quartz(width=7, height=5)
par(mfrow=c(1,3), oma=c(0,2,1,2), mar=c(5,3,3,0))
ylims <- c(0,330)
barplot (BAsame[1,], ylim=ylims, space=1, col="white")
barplot (BAdead[1,], add=T, col="black", space=1)
mtext(side=3, text="PIPO")
mtext(side=2, text=expression(paste("Mean Focal Species Basal Area (", ft^2, acre^-1, ")")), padj=-2)
barplot(BAsame[2,], ylim=ylims, yaxt="n", space=1, col="white")
barplot (BAdead[2,], add=T, col="black", space=1)
mtext(side=3, text="POTR")
barplot(BAsame[3,], ylim=ylims, yaxt="n", space=1, col="white")
barplot (BAdead[3,], add=T, col="black", space=1)
mtext(side=3, text="ABLA")

#-=--------- Log Focal species BA and BA dead ------------
# doesn't look very good, because results in some negative numbers for PIPO
# also, just generally really deceptive
quartz(width=7, height=5)
par(mfrow=c(1,3), oma=c(0,2,1,2), mar=c(5,3,3,0))
ylims <- c(0,5.5)
bas <- c(1,10,50,100,150,200)
barplot (BAln[1,], ylim=ylims, space=1, col="white", yaxt="n")
barplot (log(c(1,BAdead[1,2],1)), add=T, col="black", space=1, yaxt="n")
axis(2,at=log(bas), labels=bas)
mtext(side=3, text="PIPO")
mtext(side=2, text=expression(paste("Mean Focal Species Basal Area (", ft^2, acre^-1, ")")), padj=-2)
barplot(BAln[2,], yaxt="n", ylim=ylims, space=1, col="white", yaxt="n")
barplot (log(BAdead[2,]), add=T, col="black", space=1, yaxt="n")
axis(2,at=log(bas), labels=bas)
mtext(side=3, text="POTR")
barplot(BAln[3,], yaxt="n", ylim=ylims, space=1, col="white", yaxt="n")
barplot (log(BAdead[3,]), add=T, col="black", space=1, yaxt="n")
axis(2,at=log(bas), labels=bas)
mtext(side=3, text="ABLA")


#========================================
#          Looking at Regen 
#========================================
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

regLtmp <- by(spLreg[,columns], spLreg$Band, colMeans) # get the mean per plot regen
mspLreg <- rbind(rep(NA, time=4), regLtmp[[2]], regLtmp[[3]]) # put it into rows = L, M, Hcolumns=regen classes
regMtmp <- by(spMreg[,6:9], spMreg$Band, colMeans)
mspMreg <- rbind(regMtmp[[1]], regMtmp[[2]], regMtmp[[3]])
regHtmp <- by(spHreg[,6:9], spHreg$Band, colMeans, na.rm=TRUE)
mspHreg <- rbind(regHtmp[[1]], regHtmp[[2]], regHtmp[[3]])

mspLreg <- mspLreg / 5^2*pi * 10000  # changing to saplings/ ha
mspMreg <- mspMreg / 5^2*pi * 10000  # in saplings/plot / m2/plot * m2/ha
mspHreg <- mspHreg / 5^2*pi * 10000 


require(RColorBrewer)
quartz(width=7, height=5)
par(mfrow=c(1,3), oma=c(0,2,1,2), mar=c(5,3,3,0))
# plotting specifications
cols = brewer.pal(4,"Blues")
names <- c("Low", "Mid", "High")
ylims <- c(0, 31000)
cex.axis <- 1.2
cex.names <-1.1
ylabel <- "Saplings per ha"
barplot(t(mspLreg), beside=T, ylim=ylims, col=cols, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, legend.text = c("size 1", "size 2", "size 3", "size 4"))
mtext(side=2, text=ylabel, padj=-2.6)
barplot(t(mspMreg), beside=T, ylim=ylims, col=cols, names.arg=names, cex.names=cex.names, cex.axis=cex.axis)
barplot(t(mspHreg), beside=T, ylim=ylims, col=cols, names.arg=names, cex.names=cex.names, cex.axis=cex.axis)


#===================================================
#           Plotting abundance from my transects
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


# with all species
quartz(width=7, height=5)
cols <- brewer.pal(n=7, name="Set1")
barplot(mAbunRel, beside=T, names.arg=colnames(mAbunRel), border=cols, legend.text=rownames(mAbunRel), col=cols)


#############################
# using a loess smoothed thing instead


### melting full abundance dataset for loess smoothing
tmp <- t(mAbunRel)  
require(reshape)
mAbunRel.df <- melt (tmp,variable_name=colnames(tmp))
names(mAbunRel.df) <- c("Elev", "Species", "RelAbun")


quartz(width=7.5, height=5)
span <- 0.7
cols <- brewer.pal(n=7, name="Set1")
cols[6] <- brewer.pal(n=8, name="Accent")[8]
cols[c(2,3)] <- cols[c(3,2)]
plot(RelAbun~Elev, data=mAbunRel.df, type="n", ylab="Mean Relative Abundance", xlab= "Elevation (m)", yaxs="i", ylim=c(0,1.1), xlim=c(150,1550), mgp = c(3,1.7,0))
for (i in 1:length(unique(mAbunRel.df$Species)))
{
  sp <- unique(mAbunRel.df$Species)[i]
  test <- loess(RelAbun~Elev, data=mAbunRel.df[mAbunRel.df$Species==sp,],span=span)
  lines(test$fitted ~ mAbunRel.df$Elev[mAbunRel.df$Species==sp], lwd=3, col=cols[i])
}
par(font=3)
#legend("topleft", legend =unique(mAbunRel.df$Species), col=cols, lwd=3, bty="n" )
#legend("topright",legend=c("P. ponderosa", "P. termuloides", "A. lasiocarpa", "Picea sp.", "P. edulis", "Juniperus sp.", "P. menziesii"),col=cols[c(5,3,1,2,7,6,4)], lwd=3, bty="n", cex=.9)
par(font=1)

xlocs <- c(210.6710,  384.9774,  629.4707, 1164.3944,  863.4570, 1379.1790, 330.6715)
ylocs <- c(0.5222512, 0.2242332, 1.0154304, 0.7749657, 0.3172240, 0.1679085, 0.098828)
text(x=xlocs, y=ylocs, labels=c("Tsuga\nheterophylla","Thuja plicata", "Pseudotsuga menziesii", "Abies\nlasiocarpa", "Abies grandis", "Abies\namabilis", "Acer\n macrophyllum"),col=cols, font=3,pos=4, cex=c(1,.7,1,1,.7,.7,.7))
text(x=c(150,1550), y=c(0.4,0.4), labels = c("Elwha River", "Tree Line"), srt=90 )

#==============================================================



##################################################
#    Putting it together in 1 figure      
##################################################
display.brewer.all() # for picking a palette

### ATTEMPT 1 ####################
# mat <- matrix(c(1,2,3,
#                 4,4,4,
#                 5,6,7),nrow=3,ncol=3, byrow=T)     #create a matrix with four rows and one column, plot in order 1,2,3,4
# mat
# quartz(width=6, height=6.9)
# par(mar=c(3,5,0,0), oma=c(0,0,1,1))
# cols <- brewer.pal(n=7, name="Set1")
# #vector of heights for each plot, vector of widths for each plot
# #heights proportional to maximum catch in each group
# layout(mat=mat,heights=c(.8,.8,.8, 1, .8,.8,.8))
# layout.show(7)  #shows what the figures borders are 
# 
# #================ Plots 1:3 = BA and Mort =================
# #quartz(width=7, height=5)
# #par(mfrow=c(1,3), mar=c(5,3,3,0))
# ylims <- c(0,170)
# barplot (BAsame[1,], ylim=ylims, space=1, col=paste(cols[5], 22, sep=""), border=cols[5])
# barplot (BAdead[1,], add=T, col="black", space=1)
# #mtext(side=3, text="PIPO")
# mtext(side=2, text=expression(paste("Mean Basal Area (", ft^2, acre^-1, ")")), padj=-2)
# barplot(BAsame[2,], ylim=ylims, yaxt="n", space=1, col=paste(cols[3], 22, sep=""), border=cols[3])
# barplot (BAdead[2,], add=T, col="black", space=1)
# #mtext(side=3, text="POTR")
# barplot(BAsame[3,], ylim=ylims, yaxt="n", space=1, col=paste(cols[1], 22, sep=""), border=cols[1])
# barplot (BAdead[3,], add=T, col="black", space=1)
# #mtext(side=3, text="ABLA")
# 
# #================ Plot 4 = Abundance =================
# 
# span <- 0.5
# plot(RelAbun~Elev, data=mAbunRel.df, type="n", ylab="Relative Abundance", xlab= "Elevation (m)", yaxs="i", ylim=c(0,1.1))
# for (i in 1:length(unique(mAbunRel.df$Species)))
# {
#   sp <- unique(mAbunRel.df$Species)[i]
#   test <- loess(RelAbun~Elev, data=mAbunRel.df[mAbunRel.df$Species==sp,],span=span)
#   lines(test$fitted ~ mAbunRel.df$Elev[mAbunRel.df$Species==sp], lwd=3, col=cols[i])
# }
# par(font=3)
# legend("topright",legend=c("P. ponderosa", "P. termuloides", "A. lasiocarpa", "Picea sp.", "P. edulis", "Juniperus sp.", "P. menziesii"),col=cols[c(5,3,1,2,7,6,4)], lwd=3, bty="n", cex=.9)
# par(font=1)
# 
# 
# #================= Plots 5:7 = regen =============
# names <- c("Low", "Mid", "High")
# ylims <- c(0, 26000)
# cex.axis <- 1.2
# cex.names <-1.1
# ylabel <- "Saplings per ha"
# cols = brewer.pal(4,"Oranges")
# barplot(t(mspLreg), beside=T, ylim=ylims, col=cols, names.arg=names, cex.names=cex.names, cex.axis=cex.axis)#, legend.text = c("size 1", "size 2", "size 3", "size 4"))
# mtext(side=2, text=ylabel, padj=-2.6)
# cols = brewer.pal(4, "Greens")
# barplot(t(mspMreg), beside=T, ylim=ylims, col=cols, names.arg=names, cex.names=cex.names, cex.axis=cex.axis)
# cols= brewer.pal(4, "Reds")
# barplot(t(mspHreg), beside=T, ylim=ylims, col=cols, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, legend.text = c("size 1", "size 2", "size 3", "size 4"))





#==========================================================
#        Attempt 2 (moving things around, spacing, colors)
#                    Final Figure for R Graphics
#==========================================================
# makes a figure with BA on top, regen in middle, and total abundance below

# code for putting transect on the bottom and mort and regen above it

quartz(width=6.5, height=6.7)
mat <- matrix(c(1,2,3,
                4,5,6,
                7,7,7),nrow=3,ncol=3, byrow=T) 
par(mar=c(4,5,0,2), oma=c(2,2,3,0))
cols <- brewer.pal(n=6, name="Set1")
cols[6] <- brewer.pal(n=8, name="Accent")[8]
#vector of heights for each plot, vector of widths for each plot
#heights proportional to maximum catch in each group
layout(mat=mat,heights=c(.8,.8,.8, .8, .8,.8,1))
layout.show(7)  #shows what the figures borders are 

cex.axis <- 1.1
cex.names <-1.1
cex.lab <- .8
#================ Plots 1:3 = BA and Mort =================
#quartz(width=7, height=5)
#par(mfrow=c(1,3), mar=c(5,3,3,0))
# is initially in ft2 per acre, need to convert to m2 per ha
# ft2/acre * 0.092903m2/ft * 1acre/0.404686 ha
names <- c("Low", "Mid", "High")
ylims <- c(0,330* 0.092903/0.404686) # convert things into m^2/ha rather than ft2/acre
barplot ((BAsame[1,] * 0.092903/0.404686), ylim=ylims, space=1, col=paste(cols[1], 22, sep=""), border=cols[1], names.arg=names, cex.names=cex.names, cex.axis=cex.axis)
barplot (BAdead[1,]* 0.092903/0.404686, add=T, col="black", space=1, border=cols[1], xaxt="n", yaxt="n")
legend("topleft", legend=c("Live", "Dead"), fill=c(paste(cols[1], 22, sep=""), "black"), border=cols[1], bty="n", cex=.9)
mtext(side=3, text="Tsuga heterophylla", font=3, padj=-0.6, cex=cex.lab)
mtext(side=2, text=expression(paste("Mean Basal Area (", m^2, ha^-1, ")")), line=3.1)
barplot(BAsame[2,]* 0.092903/0.404686, ylim=ylims, yaxt="n", space=1, col=paste(cols[2], 22, sep=""), border=cols[2], names.arg=names, cex.names=cex.names, cex.axis=cex.axis)
barplot (BAdead[2,]* 0.092903/0.404686, add=T, col="black", space=1, border=cols[2], xaxt="n")
mtext(side=3, text="Pseudotsuga menziesii", font=3, padj=-0.6, cex=cex.lab)
barplot(BAsame[3,]* 0.092903/0.404686, ylim=ylims, yaxt="n", space=1, col=paste(cols[3], 22, sep=""), border=cols[3], names.arg=names, cex.names=cex.names, cex.axis=cex.axis)
barplot (BAdead[3,]* 0.092903/0.404686, add=T, col="black", space=1, border=cols[3], xaxt="n")
mtext(side=3, text="Abies lasiocarpa", font=3, padj=-0.6, cex=cex.lab)


#================= Plots 4:6 = regen =============
names <- c("Low", "Mid", "High")
ylims <- c(0, 15000)
ylabel <- "Saplings per ha"
colstmp = brewer.pal(4,"Reds")
barplot(t(mspLreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[1])#, legend.text = c("size 1", "size 2", "size 3", "size 4"))
mtext(side=2, text=ylabel,line=3.1)
legend("topleft", legend=c("<0.5m height", "<1.3m height", "<2.5cm DBH", "2.5-5cm DBH"), fill=colstmp, border=cols[1], bty="n", cex=.8)
colstmp = brewer.pal(4, "Blues")
barplot(t(mspMreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[2])
colstmp= brewer.pal(4, "Greens")
barplot(t(mspHreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[3])

#================ Plot 7 = Abundance =================

span <- 0.7
cols <- brewer.pal(n=7, name="Set1")
cols[6] <- brewer.pal(n=8, name="Accent")[8]

cols[c(2,3,4)] <- cols[c(4,2,3)]
plot(RelAbun~Elev, data=mAbunRel.df, type="n", ylab="Mean Relative Abundance", xlab= "Elevation (m)", yaxs="i", ylim=c(0,1.1), xlim=c(150,1550), mgp = c(3,1.7,0))
for (i in 1:length(unique(mAbunRel.df$Species)))
{
  sp <- unique(mAbunRel.df$Species)[i]
  test <- loess(RelAbun~Elev, data=mAbunRel.df[mAbunRel.df$Species==sp,],span=span)
  lines(test$fitted ~ mAbunRel.df$Elev[mAbunRel.df$Species==sp], lwd=3, col=cols[i])
}
par(font=3)
#legend("topleft", legend =unique(mAbunRel.df$Species), col=cols, lwd=3, bty="n" )
#legend("topright",legend=c("P. ponderosa", "P. termuloides", "A. lasiocarpa", "Picea sp.", "P. edulis", "Juniperus sp.", "P. menziesii"),col=cols[c(5,3,1,2,7,6,4)], lwd=3, bty="n", cex=.9)
par(font=1)

xlocs <- c(210.6710,  384.9774,  629.4707, 1164.3944,  863.4570, 1379.1790, 330.6715)
ylocs <- c(0.5222512, 0.2242332, 1.0154304, 0.7749657, 0.3172240, 0.1679085, 0.098828)
text(x=xlocs, y=ylocs, labels=c("Tsuga\nheterophylla","Thuja plicata", "Pseudotsuga menziesii", "Abies\nlasiocarpa", "Abies grandis", "Abies\namabilis", "Acer\n macrophyllum"),col=cols, font=3,pos=4, cex=c(1,.7,1,1,.7,.7,.7))
text(x=c(150,1550), y=c(0.4,0.4), labels = c("Elwha River", "Tree Line"), srt=90 )





span <- 0.6
cols <- brewer.pal(n=7, name="Set1")
cols[6] <- brewer.pal(n=8, name="Accent")[8]
cols[7] <- brewer.pal(n=8, name="Set1")[8]
cols <- cols[c(1,4,5,6,2,3,7)]
plot(RelAbun~Elev, data=mAbunRel.df, type="n", ylab="Relative Abundance", xlab= "Elevation (m)", yaxs="i", ylim=c(0,1.1), xlim=c(950, 2250), mgp = c(3,1.7,0))
for (i in 1:length(unique(mAbunRel.df$Species)))
{
  sp <- unique(mAbunRel.df$Species)[i]
  test <- loess(RelAbun~Elev, data=mAbunRel.df[mAbunRel.df$Species==sp,],span=span, degree=2)
  lines(test$fitted ~ mAbunRel.df$Elev[mAbunRel.df$Species==sp], lwd=3, col=cols[i])
}
par(font=3)
#legend("topright", legend =unique(mAbunRel.df$Species), col=cols, lwd=3, bty="n" )
legend(x=1870, y=0.4,legend=c("Thuja plicata", "Larix occidentalis", "Pinus spp.", "Picea engelmannii"),col=cols[c(2,3,4,7)], lwd=3, bty="n", cex=.8)
par(font=1)
#locs <- locator(n=3)
locs$x[2] <- 1270
text(x=locs$x, y=locs$y, labels=c("Tsuga\nheterophylla", "Pseudotsuga\nmenziesii", "Abies\nlasiocarpa"),col=cols[c(1,5,6)], font=3,pos=4)
text(x=c(950,2250), y=c(0.4,0.4), labels = c("Lake McDonald", "Tree Line"), srt=90 )


