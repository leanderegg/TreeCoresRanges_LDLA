###################################################
#       Plotting Relative Abundance for NP reporting
#                For WA transect
###################################################
# LDL Anderegg
# 3/29/14
# code lifted directly from bottom of CO_assumptionchecks.R


require(RColorBrewer)
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







