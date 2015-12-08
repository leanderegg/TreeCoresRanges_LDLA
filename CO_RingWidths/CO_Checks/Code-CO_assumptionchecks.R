#####################################
####### Testing CO assumptions ######
############# 12/23/13 ##############


# read in scratch csv for looking at PIPO
COtrees <- read.table("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_Checks/CO_prelim_aspectDBH_12-23-13.csv", head=T, sep=",")
COtreesnew <- read.csv("/Users/leeanderegg/Dropbox/treecoresranges(leander) (1)/CO_RingWidths/CO_Checks/CO_focaltree_3_6.csv", header=T)
      # good. all trees seem to be there.
# checking the aspect
aspect <- table(COtrees$Aspect,COtrees$Band,COtrees$Species)
PIPOasp <- data.frame(aspect[,,"PIPO"])
POTRasp <- data.frame(aspect[,,"POTR"])
ABLAasp <- data.frame(aspect[,,"ABLA"])
degrees <- c(180,135,157.5,202.5,225,270,247.2)

## PIPO
PH <- rep(x=degrees, times= PIPOasp$H) # repeating the value of degrees the number of times dictated by the column from my table
PM <- rep  (x=degrees, times= PIPOasp$M)
PL <- rep (x=degrees, times= PIPOasp$L)
(masp<- c(mean(PH), mean(PM), mean(PL))) # seems like everything is ok here
t.test(x=PM,y=PL)
par(mfrow=c(3,1))
hist(PH, ylab="High trees")
hist(PM, ylab="Med trees")
hist(PL, ylab="Low trees")


## POTR
POH <- rep(x=degrees, times= POTRasp$H) # repeating the value of degrees the number of times dictated by the column from my table
POM <- rep  (x=degrees, times= POTRasp$M)
POL <- rep (x=degrees, times= POTRasp$L)
(masp<- c(mean(POH), mean(POM), mean(POL))) # seems like everything is ok here
hist(POH, ylab="High trees")
hist(POM, ylab="Med trees")
hist(POL, ylab="Low trees")

## ABLA
AH <- rep(x=degrees, times= ABLAasp$H) # repeating the value of degrees the number of times dictated by the column from my table
AM <- rep  (x=degrees, times= ABLAasp$M)
AL <- rep (x=degrees, times= ABLAasp$L)
(masp<- c(mean(AH), mean(AM), mean(AL))) 
quartz(width=4, heigh=8)
par(mfrow=c(3,1))
hist(AH)
hist(AM)
hist(AL) 
#meh, the means seem to be different,
  #but the histograms look like I did an ok job staying within the same aspect



######## DBH across elevation ############
PIDBH <- COtrees[COtrees$Species=="PIPO",]
H <- PIDBH$DBH[PIDBH$Band=="H"]  
M <- PIDBH$DBH[PIDBH$Band=="M"]
L <- PIDBH$DBH[PIDBH$Band=="L"]
(mDBH <- c(mean(PIDBH$DBH[PIDBH$Band=="H"]), mean(PIDBH$DBH[PIDBH$Band=="M"]), mean(PIDBH$DBH[PIDBH$Band=="L"])))
summary(aov(PIDBH$DBH~PIDBH$Band))
boxplot(PIDBH$DBH~PIDBH$Band)
# well, high elevation has signifcantly larger trees than low 
  # (maybe due to the logging, but also possibly because of better growth?)

PODBH <- COtrees[COtrees$Species=="POTR",]
boxplot(DBH~Band, data=PODBH)
summary(aov(DBH~Band, data= PODBH))
# again, looks like high elevation has significantly larger trees than low
  # this one has to be attributed to smaller trees at low due to physiological constraints

ABDBH <- COtrees[COtrees$Species=="ABLA",]
boxplot(DBH~Band, ABDBH)
# low trees are bigger than high here, which is totally different
  # from what I would expect with high trees being unlogged and low and mid logged.



#=========== Using new, clean, and complete Dataz ===================
boxplot(DBH~Band+Species, data=COtreesnew)
COtreesnew$Height <- COtreesnew$dist *(COtreesnew$p_top.of.canopy.-COtreesnew$p_base.ref.)
boxplot(Height ~ Band + Species, data=COtreesnew)
# supports what I was thinking for the DBH trends above
boxplot(BA_tot~Band+Species+Comp, data=COtreesnew)

levels(COtreesnew$Band) <- list(L = "L", M="M", H="H") # reordering the levels so they plot low to high
levels(COtreesnew$Species) <- list(PIPO="PIPO", POTR="POTR", ABLA="ABLA")

# Stand Basal Area
quartz(width=7, height=5)
par(mfrow=c(1,3),mar=c(5,5,3,0), oma=c(0,1,2,1))
for(i in 1:3){
  sp <- levels(COtreesnew$Species)[i]
  boxplot(BA_tot~Comp + Band,at=c(1,2,4,5,7,8), xaxt="n", data=COtreesnew[which(COtreesnew$Species==sp),], main= sp, col=c("grey", "white"), ylim=c(0,350) )#, ylab=expression(paste("Stand Basal Area (", ft^2, ha^-1, ")")))
  axis(side=1,at=c(1.5,4.5,7.5), labels=c("Low", "Mid", "High"), cex.axis=1.4)
}
mtext(text=expression(paste("Stand Basal Area (", ft^2, ha^-1, ")")),side=2,outer=T, line=-2, cex=1.5)


# DBH
quartz(width=7, height=5)
par(mfrow=c(1,3),mar=c(5,5,3,0), oma=c(0,1,2,1))
for(i in 1:3){
  sp <- levels(COtreesnew$Species)[i]
  boxplot(DBH~Comp + Band,at=c(1,2,4,5,7,8), xaxt="n", data=COtreesnew[which(COtreesnew$Species==sp),], main= sp, col=c("grey", "white"), ylim=c(0, 80) )#, ylab=expression(paste("Stand Basal Area (", ft^2, ha^-1, ")")))
  axis(side=1,at=c(1.5,4.5,7.5), labels=c("Low", "Mid", "High"), cex.axis=1.4)
}
mtext(text="Tree DBH",side=2,outer=T, line=-2, cex=1.5)


# slope
quartz(width=7, height=5)
par(mfrow=c(1,3),mar=c(5,5,3,0), oma=c(0,1,2,1))
for(i in 1:3){
  sp <- levels(COtreesnew$Species)[i]
  boxplot(Slope~Comp + Band,at=c(1,2,4,5,7,8), xaxt="n", data=COtreesnew[which(COtreesnew$Species==sp),], main= sp, col=c("grey", "white") , ylim=c(0, 50) )#, ylab=expression(paste("Stand Basal Area (", ft^2, ha^-1, ")")))
  axis(side=1,at=c(1.5,4.5,7.5), labels=c("Low", "Mid", "High"), cex.axis=1.4)
}
mtext(text="Slope (%)",side=2,outer=T, line=-2, cex=1.5)


# Conspecific live BA
quartz(width=7, height=5)
par(mfrow=c(1,3),mar=c(5,5,3,0), oma=c(0,1,2,1))
for(i in 1:3){
  sp <- levels(COtreesnew$Species)[i]
  boxplot(BA_sameL~Comp + Band,at=c(1,2,4,5,7,8), xaxt="n", data=COtreesnew[which(COtreesnew$Species==sp),], main= sp, col=c("grey", "white") , ylim=c(0, 350) )#, ylab=expression(paste("Stand Basal Area (", ft^2, ha^-1, ")")))
  axis(side=1,at=c(1.5,4.5,7.5), labels=c("Low", "Mid", "High"), cex.axis=1.4)
}
mtext(text="Conspecific live BA",side=2,outer=T, line=-2, cex=1.5)

# Trees within 5m
quartz(width=7, height=5)
par(mfrow=c(1,3),mar=c(5,5,3,0), oma=c(0,1,2,1))
for(i in 1:3){
  sp <- levels(COtreesnew$Species)[i]
  boxplot(in5_tot~Comp + Band,at=c(1,2,4,5,7,8), xaxt="n", data=COtreesnew[which(COtreesnew$Species==sp),], main= sp, col=c("grey", "white") , ylim=c(0, 28) )#, ylab=expression(paste("Stand Basal Area (", ft^2, ha^-1, ")")))
  axis(side=1,at=c(1.5,4.5,7.5), labels=c("Low", "Mid", "High"), cex.axis=1.4)
}
mtext(text="Trees in 5m",side=2,outer=T, line=-2, cex=1.5)


# Number of Crowns
quartz(width=7, height=5)
par(mfrow=c(1,3),mar=c(5,5,3,0), oma=c(0,1,2,1))
for(i in 1:3){
  sp <- levels(COtreesnew$Species)[i]
  boxplot(N_Cr~Comp + Band,at=c(1,2,4,5,7,8), xaxt="n", data=COtreesnew[which(COtreesnew$Species==sp),], main= sp, col=c("grey", "white") , ylim=c(0, 8) )#, ylab=expression(paste("Stand Basal Area (", ft^2, ha^-1, ")")))
  axis(side=1,at=c(1.5,4.5,7.5), labels=c("Low", "Mid", "High"), cex.axis=1.4)
}
mtext(text="Number of touching crowns",side=2,outer=T, line=-2, cex=1.5)


# Height
COtreesnew$Height <- (COtreesnew$p_top.of.canopy. - COtreesnew$p_base.ref.) * COtreesnew$dist/100
quartz(width=7, height=5)
par(mfrow=c(1,3),mar=c(5,5,3,0), oma=c(0,1,2,1))
for(i in 1:3){
  sp <- levels(COtreesnew$Species)[i]
  boxplot(Height~Comp + Band,at=c(1,2,4,5,7,8), xaxt="n", data=COtreesnew[which(COtreesnew$Species==sp),], main= sp, col=c("grey", "white") , ylim=c(5, 30) )#, ylab=expression(paste("Stand Basal Area (", ft^2, ha^-1, ")")))
  axis(side=1,at=c(1.5,4.5,7.5), labels=c("Low", "Mid", "High"), cex.axis=1.4)
}
mtext(text="Tree Height (m)",side=2,outer=T, line=-2, cex=1.5)




bandtest <- COtreesnew$Band
levels(bandtest)<- list(L = "L", M="M", H="H")
# ok, so high and low competition comparisons are not quite as straight forward as I would like,
# because the high and low treatments pretty much vary between elevations for all species

#===============================================
#====== checking out mortality roughly =========
quartz(width=7, height=5)
par(mfrow=c(1,3),mar=c(5,5,3,0))
for(i in 1:3){
  sp <- levels(COtreesnew$Species)[i]
  boxplot(BA_totR~ Band, data=COtreesnew[which(COtreesnew$Species==sp),], 
          main= sp, col=c("grey", "white"), 
          ylab=expression(paste("Stand Basal Area (", ft^2, ha^-1, ")")))
}

## Making matrix of average % BA mortality for each species & Band
mort <- COtreesnew$BA_sameRG/COtreesnew$BA_same_all # getting percentage dead of total same basal area
mort[which(is.na(mort))] <- 0
perc_mort <- tapply(mort, INDEX=list(COtreesnew$Species, COtreesnew$Band), FUN=mean)
perc_mort <- perc_mort[,c(2,3,1)] # reorder elevations to go Low Med High
perc_mort <- perc_mort[c(2,3,1),] # reorder species to go PIPO, POTR, ABLA

#plotting
quartz(width=7, height=5)
par(mfrow=c(1,3), oma=c(0,2,1,2), mar=c(5,3,3,0))
barplot (perc_mort[1,], ylim=c(0,0.30), space=1)
mtext(side=3, text="PIPO")
mtext(side=2, text= "% basal area dead", padj=-2.7)
barplot(perc_mort[2,], ylim=c(0, 0.3), yaxt="n", space=1)
mtext(side=3, text="POTR")
barplot(perc_mort[3,], ylim=c(0, 0.3), yaxt="n", space=1)
mtext(side=3, text="ABLA")


#### Making barplot of mean BA with mortality colored in ###
BAsame <- tapply(COtreesnew$BA_same_all, INDEX=list(COtreesnew$Species, COtreesnew$Band), FUN=mean)
BAsame <- BAsame[,c(2,3,1)] # reorder elevations to go Low Med High
BAsame <- BAsame[c(2,3,1),]
BAln <- log(BAsame)

BAdead <- BAsame * perc_mort # and make a vector that's the average dead BA

#-=--------- Focal species BA and BA dead ------------
quartz(width=7, height=5)
par(mfrow=c(1,3), oma=c(0,2,1,2), mar=c(5,3,3,0))
ylims <- c(0,170)
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


require(RColorBrewer)
quartz(width=7, height=5)
par(mfrow=c(1,3), oma=c(0,2,1,2), mar=c(5,3,3,0))
# plotting specifications
cols = brewer.pal(4,"Blues")
names <- c("Low", "Mid", "High")
ylims <- c(0, 26000)
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
abun <- read.csv(file="/Users/leeanderegg/Dropbox/treecoresranges(leander) (1)/CO_RingWidths/CO_Checks/CO_abundancetransects3_7.csv", header=T)

#limiting down to only live trees
abunlive <- abun[,-grep("d", x=names(abun), ignore.case=FALSE)]

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

write.csv(mAbunRel, "CO_Relative_Abundance.csv")

# with all species
quartz(width=7, height=5)
cols <- brewer.pal(n=7, name="Set1")
barplot(mAbunRel, beside=T, names.arg=colnames(mAbunRel), border=cols, legend.text=rownames(mAbunRel), col=cols)

# combining PSME and Junipers
other <- colSums(mAbunRel[c("PSME","Juniperus_sp"),])
mAbun2 <- rbind(mAbunRel[c(1,2,3,5,7),], other)
quartz(width=7, height=5)
cols <- brewer.pal(n=6, name="Set1")
barplot(mAbun2, beside=T, names.arg=colnames(mAbun2), border=cols, legend.text=rownames(mAbun2), col=cols)
#barplot(mAbun2, beside=F, names.arg=colnames(mAbun2), border=cols, legend.text=rownames(mAbun2), col=cols)

#############################
# using a loess smoothed thing instead


### melting full abundance dataset for loess smoothing
tmp <- t(mAbunRel)  
require(reshape)
mAbunRel.df <- melt (tmp,variable_name=colnames(tmp))
names(mAbunRel.df) <- c("Elev", "Species", "RelAbun")


quartz(width=7, height=5)
span <- 0.5
cols <- cols <- brewer.pal(n=7, name="Set1")
plot(RelAbun~Elev, data=mAbunRel.df, type="n", ylab="Relative Abundance", xlab= "Elevation (m)", yaxs="i", ylim=c(0,1.1))
for (i in 1:length(unique(mAbunRel.df$Species)))
{
  sp <- unique(mAbunRel.df$Species)[i]
  test <- loess(RelAbun~Elev, data=mAbunRel.df[mAbunRel.df$Species==sp,],span=span)
  lines(test$fitted ~ mAbunRel.df$Elev[mAbunRel.df$Species==sp], lwd=3, col=cols[i])
}
par(font=3)
legend("topright",legend=c("P. ponderosa", "P. termuloides", "A. lasiocarpa", "Picea sp.", "P. edulis", "Juniperus sp.", "P. menziesii"),col=cols[c(5,3,1,2,7,6,4)], lwd=3, bty="n", cex=.9)
par(font=1)


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

# first getting rid of PSME, because you can't see it
mAbunRel.df2 <- mAbunRel.df[which(mAbunRel.df$Species!="PSME"),]
levels(mAbunRel.df$Species) <- list(PIPO = "PIPO", POTR="POTR", ABLA="ABLE", PIEN = "Picea_sp", PIED = "PIED", JUNI = "Juniperus_sp")


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

