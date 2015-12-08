###################################################
#       Plotting Relative abundance for NP reporting
#                For MT transect
###################################################
# LDL Anderegg
# 3/29/14
# code lifted directly from bottom of CO_assumptionchecks.R


### Things to note: I only measured abundance transects up to 1700m, though my transect went up to 2260 ish
# so to make the abundance figure, I had to continue graph out based on either BAI data from ABLA M and H or w/in 5m from the same

require(RColorBrewer)
#===================================================
#           Plotting abundance from my transects
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


# with all species
#quartz(width=7, height=5)
#cols <- brewer.pal(n=8, name="Set1")
#barplot(mAbunRel, beside=T, names.arg=colnames(mAbunRel), border=cols, legend.text=rownames(mAbunRel), col=cols)


#############################
# using a loess smoothed thing instead


### melting full abundance dataset for loess smoothing
tmp <- t(mAbunRel)  
require(reshape)
mAbunRel.df <- melt (tmp,variable_name=colnames(tmp))
names(mAbunRel.df) <- c("Elev", "Species", "RelAbun")


quartz(width=7.5, height=5)
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
legend(x=1850, y=0.3,legend=c("Thuja plicata", "Larix occidentalis", "Pinus spp.", "Picea engelmannii"),col=cols[c(2,3,4,7)], lwd=3, bty="n", cex=.8)
par(font=1)
#locs <- locator(n=3)
locs$x[2] <- 1270
text(x=locs$x, y=locs$y, labels=c("Tsuga\nheterophylla", "Pseudotsuga\nmenziesii", "Abies\nlasiocarpa"),col=cols[c(1,5,6)], font=3,pos=4)
text(x=c(950,2250), y=c(0.4,0.4), labels = c("Lake McDonald", "Tree Line"), srt=90 )

