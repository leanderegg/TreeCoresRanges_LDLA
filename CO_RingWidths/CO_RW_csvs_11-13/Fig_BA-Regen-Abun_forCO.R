

#==========================================================
#        Figure 1 (Tree BA, Regeneration, and Relative Abundance)
#                    Final Figure for R Graphics
#==========================================================
# makes a figure with BA on top, regen in middle, and total abundance below

require(RColorBrewer)

#======================= Reading in Data =============================

###### Reading in Relative Abundance Data
mAbunRel <- read.csv("CO_Relative_Abundance.csv", header=T)

# melting full abundance dataset for loess smoothing
tmp <- t(mAbunRel)  
require(reshape)
mAbunRel.df <- melt (tmp,variable_name=colnames(tmp))
names(mAbunRel.df) <- c("Elev", "Species", "RelAbun")

# first getting rid of PSME, because you can't see it
mAbunRel.df2 <- mAbunRel.df[which(mAbunRel.df$Species!="PSME"),]
levels(mAbunRel.df$Species) <- list(PIPO = "PIPO", POTR="POTR", ABLA="ABLE", PIEN = "Picea_sp", PIED = "PIED", JUNI = "Juniperus_sp")



###### Reading in Basal Area data
# full focal tree database
COtreesnew <- read.csv("CO_focaltree_3_6.csv", header=T)

# pulling out Basal Area of the same species
BAsame <- tapply(COtreesnew$BA_same_all, INDEX=list(COtreesnew$Species, COtreesnew$Band), FUN=mean)
BAsame <- BAsame[,c(2,3,1)] # reorder elevations properly
BAsame <- BAsame[c(2,3,1),] # reorder species to go PIPO, POTR, ABLA
BAdead <- BAsame * perc_mort # and make a vector that's the average dead BA


## Making matrix of average % BA mortality for each species & Band
mort <- COtreesnew$BA_sameRG/COtreesnew$BA_same_all # getting percentage dead of total same basal area
mort[which(is.na(mort))] <- 0
perc_mort <- tapply(mort, INDEX=list(COtreesnew$Species, COtreesnew$Band), FUN=mean)
perc_mort <- perc_mort[,c(2,3,1)] # reorder elevations to go Low Med High
perc_mort <- perc_mort[c(2,3,1),] # reorder species to go PIPO, POTR, ABLA


######## Formatting and converting Regeneration Data
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

mspLreg <- mspLreg / 5^2*pi * 10000  # changing to saplings/ ha
mspMreg <- mspMreg / 5^2*pi * 10000  # in saplings/plot / m2/plot * m2/ha
mspHreg <- mspHreg / 5^2*pi * 10000 



#================== Plot Code ===================================

# Makes a multipanel plot with BA barplots on top, Regen barplots in the middle
# and one large panel of relative Abun on the bottom


quartz(width=6.5, height=6.7)
# set up the matrix for my layout
mat <- matrix(c(1,2,3,
                4,5,6,
                7,7,7),nrow=3,ncol=3, byrow=T) 

# split up my plotting device as specified above
layout(mat=mat,heights=c(.8,.8,.8, .8, .8,.8,1))
layout.show(7)  #shows what the figures borders are 

#set my par preferences
par(mar=c(4,5,0,2), oma=c(2,2,3,0))

#set some colors from different pallettes
cols <- brewer.pal(n=6, name="Set1")
cols[6] <- brewer.pal(n=8, name="Accent")[8]


cex.axis <- 1.1
cex.names <-1.1
cex.lab <- .8
#================ Plots 1:3 = BA and Mort =================
#quartz(width=7, height=5)
#par(mfrow=c(1,3), mar=c(5,3,3,0))
# is initially in ft2 per acre, need to convert to m2 per ha
# ft2/acre * 0.092903m2/ft * 1acre/0.404686 ha
names <- c("Low", "Mid", "High")
ylims <- c(0,175* 0.092903/0.404686)
barplot ((BAsame[1,] * 0.092903/0.404686), ylim=ylims, space=1, col=paste(cols[1], 22, sep=""), border=cols[1], names.arg=names, cex.names=cex.names, cex.axis=cex.axis)
barplot (BAdead[1,]* 0.092903/0.404686, add=T, col="black", space=1, border=cols[1], xaxt="n", yaxt="n")
mtext(side=3, text="Pinus ponderosa", font=3, padj=-0.6, cex=cex.lab)
mtext(side=2, text=expression(paste("Mean Basal Area (", m^2, ha^-1, ")")), line=3.1)
barplot(BAsame[2,]* 0.092903/0.404686, ylim=ylims, yaxt="n", space=1, col=paste(cols[2], 22, sep=""), border=cols[2], names.arg=names, cex.names=cex.names, cex.axis=cex.axis)
barplot (BAdead[2,]* 0.092903/0.404686, add=T, col="black", space=1, border=cols[2], xaxt="n")
mtext(side=3, text="Populus tremuloides", font=3, padj=-0.6, cex=cex.lab)
barplot(BAsame[3,]* 0.092903/0.404686, ylim=ylims, yaxt="n", space=1, col=paste(cols[3], 22, sep=""), border=cols[3], names.arg=names, cex.names=cex.names, cex.axis=cex.axis)
barplot (BAdead[3,]* 0.092903/0.404686, add=T, col="black", space=1, border=cols[3], xaxt="n")
mtext(side=3, text="Abies lasiocarpa", font=3, padj=-0.6, cex=cex.lab)
legend("topright", legend=c("Live", "Dead"), fill=c(paste(cols[3], 22, sep=""), "black"), border=cols[3], bty="n", cex=1.2, names.arg=names, cex.names=cex.names)


#================= Plots 4:6 = regen =============
names <- c("Low", "Mid", "High")
ylims <- c(0, 26000)
ylabel <- "Saplings per ha"
colstmp = brewer.pal(4,"Reds")
barplot(t(mspLreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[1])#, legend.text = c("size 1", "size 2", "size 3", "size 4"))
mtext(side=2, text=ylabel,line=3.1)
colstmp = brewer.pal(4, "Blues")
barplot(t(mspMreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[2])
colstmp= brewer.pal(4, "Greens")
barplot(t(mspHreg), beside=T, ylim=ylims, col=colstmp, names.arg=names, cex.names=cex.names, cex.axis=cex.axis, border=cols[3])
legend("topright", legend=c("<0.5m height", "<1.3m height", "<2.5cm DBH", "2.5-5cm DBH"), fill=colstmp, border=cols[3], bty="n")

#================ Plot 7 = Abundance =================

span <- 0.5
plot(RelAbun~Elev, data=mAbunRel.df2, type="n", ylab="Relative Abundance", xlab= "Elevation (m)", yaxs="i", ylim=c(0,1.1), frame.plot="n", cex.axis=cex.axis)
# can't actually get the species in the right order because of unique. need to change this up
mtext(side=2, text="Relative Abundance", line=3.1)
species <- unique(mAbunRel.df2$Species)[c(4,3,1,2,5,6)]
for (i in 1:length(species))
{
  sp <- species[i]
  test <- loess(RelAbun~Elev, data=mAbunRel.df2[mAbunRel.df2$Species==sp,],span=span)
  lines(test$fitted ~ mAbunRel.df2$Elev[mAbunRel.df2$Species==sp], lwd=3, col=cols[i])
}
par(font=3)
#legend("topright",legend=c("P. ponderosa", "P. tremuloides", "A. lasiocarpa", "Picea sp.", "P. menziesii", "Juniperus sp.", "P. edulis"),col=cols, lwd=3, bty="n", cex=.9,ncol=2)
par(font=1)
text(x=c(2410,2810, 3300,3030, 2450, 2380), y=c(0.7,0.7, 0.3, 0.3, 0.12,0.28), labels=c("P.ponderosa","P. tremuloides", "A. lasiocarpa", "Picea sp", "P. edulis", "J. osteosperma"),col=cols, font=3, cex=)
mtext(side=1, text="Elevation (m)", line=3)

