##################################################################
#######            Code for creating sparkplots   
######                for R Graphics class
######                 by: LDL Anderegg
######                 3/20/2014
##################################################################
## Requires detrended ##rwi.clim.allelev## data frame from Mixed_effects_try1.R
# as well as ##PIPO.master## from the same file

rw.clim.allelev <- read.csv("CO-PIPO-allRWIClim.csv")
PIPO.master <- read.csv("MASTER_CO-PIPO_02072014.csv")



#######################################################################
############### Final Figure for R Graphics: ##########################
############## Sparkplots plust Master chronologies ###################
#######################################################################
# puts master chronologies with individual chronologies in gray in a top plot,
# and then plots sparkplots from 20 trees (dropped 11 and 12 for high elevation)
# underneath

# set up layout
mat <- matrix(c(1,1,1,
                2,3,4),nrow=2,ncol=3, byrow=T)

quartz(width=8, height=6)
par(mar=c(3,3,1,2), oma=c(2,3,1,1))
layout(mat=mat, heights=c(.6,1,1,1))
layout.show(4)
linecol <- "#AAAAAA"
plot(PIPO.L~X, data=PIPO.master, type="n", lwd=2, col="red3", ylab="", xlab="", ylim=c(-.2, 4.7), xlim=c(1900,2013), yaxt="n", mgp = c(2,.7,0), cex.axis=1.2)
mtext(text="Detrended Master \n Chronologies",side=2, line=1.4)
text(x=c(2011,2011,2011), y=c(1.15,2.65,4.15), pos=4, labels=c("Low", "Mid", "High"), font=2)
abline(h=1, lty=2, col="gray")
abline(h=2.5, lty=2, col="gray")
abline(h=4, lty=2, col="gray")

#=========== adding individual chronologies under masters ===========
# ---------- Low Elevation gray lines --------------
toplot <- subset(x=rw.clim.allelev,subset=Elev=="L")
years <- toplot$Year
names <- toplot$Tree
data <- toplot$value
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

#----------- Mid elevation gray lines
toplot <- subset(x=rw.clim.allelev,subset=Elev=="M")
years <- toplot$Year
names <- toplot$Tree
data <- toplot$value
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
#------------ High elevation gray lines -------------
toplot <- subset(x=rw.clim.allelev,subset=Elev=="H")
years <- toplot$Year
names <- toplot$Tree
data <- toplot$value
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
#------ Add in master chronologies on top of every thing -----------
points(PIPO.L~X, data=PIPO.master, type="l", lwd=2, col="red3")#, ylab="Detrended Master Chronologies", xlab="Year", ylim=c(0, 4.5), yaxt="n")
points(PIPO.M+1.5~X, data=PIPO.master, type="l", lwd=2, col="purple3")
points(PIPO.H+3~X, data=PIPO.master, type="l", lwd=2, col="blue3")




# ================ now making sparkplots below as figs 2:4 ================
# ---------- Low Elevation Sparkplot
par(mar=c(1,0,2,1))
toplot <- subset(x=rw.clim.allelev,subset=Elev=="L")
main <- "Low elevation"
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
}
mtext(text=main, side=3, line=.4, cex=1, adj=.55)


#------------------ Mid Elevation Sparkplot -------------
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
mtext(text=main, side=3, line=.4, cex=1, adj=.55)

#-------------------- High Elevation ------------
toplot <- subset(x=rw.clim.allelev,subset=Elev=="H")
# can also knock of 11 and 12 to make each elevation have the same number of trees
toplot <- toplot[which(toplot$Tree != "H.11H" & toplot$Tree != "H.11L" & toplot$Tree != "H.12H" & toplot$Tree != "H.12L"),]
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
mtext(text=main, side=3, line=.4, cex=1, adj=.55)
