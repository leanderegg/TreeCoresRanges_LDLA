#########################################################
#########          PIPO vs Climate 
#########          Code for various figures
#########       PIPO RWI vs. climate variables
#########################################################

# Includes code for simple and more complex figures showing the relationship
# between Ponderosa pine RWI collected by LDLA in the summer of 2013 and various
# climate variables acquired from the climateWNA downscaling of the PRISM dataset

# requires data from Mixed_effects_try1.R, mostly rw.clim.allelev
require(lme4)

# original data exploration using lattice
require(lattice)
#xyplot(value~CMDGS | Elev, groups=Comp, data=rw.clim.allelev, pch=19, type =c("p","r"),main=paste("PIPO-CMDGS"),lwd=3, drop.unused.levels=FALSE, col=c("red3", "darkblue"))
#xyplot(value~TminAn | Elev, groups=Comp, data=rw.clim.allelev, pch=19, type =c("p","r"),main=paste("PIPO-CMDGS"),lwd=3, drop.unused.levels=FALSE, col=c("red3", "darkblue"))
#xyplot(value~TmaxAn | Elev, groups=Comp, data=rw.clim.allelev, pch=19, type =c("p","r"),main=paste("PIPO-CMDGS"),lwd=3, drop.unused.levels=FALSE, col=c("red3", "darkblue"))
#xyplot(value~TaveAn | Elev, groups=Comp, data=rw.clim.allelev, pch=19, type =c("p","r"),main=paste("PIPO-CMDGS"),lwd=3, drop.unused.levels=FALSE, col=c("red3", "darkblue"))
#quartz(width=6, height=5)
#xyplot(value~TminAn + TaveAn + PPTGS | Elev, groups=Comp, data=rw.clim.allelev, pch=19, type =c("p","r"),main=paste("PIPO-CMDGS"),lwd=3, drop.unused.levels=FALSE, col=c(rgb(.8,0,0,0.2), rgb(0,0,.6,0.3)))

rw.clim.allraw <- read.csv(file="CO-PIPO-allRWI_rawClim.csv")
##========= Final Figure: RWI vs Tmax, CMD, Tmin =============
#                And all elevation overplotted (with legend on top)
#               (3 climate attributes in 3 panels)
#===============================================================
quartz(width=5, height=4)
a<-"09" # choose an alpha value for opacity
ylim= c(0,2.2) # delete for PIPO figure
par(mfrow=c(1,3), oma=c(1,5,1,2), mar=c(5,0,3,0))
#xyplot(value~TminAn + TaveAn + PPTGS | Elev, groups=Comp, data=rw.clim.allraw, pch=19, type =c("p","r"),main=paste("PIPO-CMDGS"),lwd=3, drop.unused.levels=FALSE, col=c(rgb(.8,0,0,0.2), rgb(0,0,.6,0.3)))
plot(value~TmaxAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], col=paste("#000099",a,sep=""), xlim=c(22,28.3), ylim=ylim, pch=19, cex=.5, ylab="Ring Width Index", xlab="Max Temp (°C)", cex.lab=1.4, font.lab=2, cex.axis=1.3)
points(value~TmaxAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),], col=paste("#990099",a,sep=""), pch=19, cex=.5)
points(value~TmaxAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),], col=paste("#990000",a,sep=""), pch=19, cex=.5)
H <- lm(value~TmaxAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),])
M <- lm(value~TmaxAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),])
L <- lm(value~TmaxAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),])
abline(H, col="#000099", lwd=3)
abline(M, col="#990099", lwd=3)
abline(L, col="#CC0000", lwd=3)
mtext(text="Ring Width Index", 2, outer=T, padj=-3, font=2, adj=.6)
mtext(text="a)", 3, adj=0.05, font=2, padj=1.5)

plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], col=paste("#000099",a,sep=""), xlim=c(220,680), ylim=ylim, pch=19, cex=.5, yaxt="n", xlab="CMD", cex.lab=1.4, font.lab=2, cex.axis=1.3)
points(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),], col=paste("#990099",a,sep=""), pch=19, cex=.5)
points(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),], col=paste("#990000",a,sep=""), pch=19, cex=.5)
H <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),])
M <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),])
L <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),])
abline(H, col="#000099", lwd=3)
abline(M, col="#990099", lwd=3)
abline(L, col="#CC0000", lwd=3)
mtext(text="b)", 3, adj=0.05, font=2, padj=1.5)

par(font=2) #can't make bold in legend function
legend(x=450, y=2.53, xpd=NA, xjust=0.5, legend=c("High Elev", "Mid Elev", "Low Elev"), lwd=3, col=c("#000099","#990099","#CC0000"),bty="n", ncol=3, cex=1.2)
# legend y= 2.33 for PIPO
par(font=1)
# x=280, y=2.0 to get it in the plot, the xjust gives me center justification on x value, default seems to be left

plot(value~TminAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], col=paste("#000099",a,sep=""), pch=19, ylim=ylim, cex=.5, yaxt="n", xlab="Min Temp (°C)", cex.lab=1.4, font.lab=2, cex.axis=1.3)
points(value~TminAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),], col=paste("#990099",a,sep=""), pch=19, cex=.5)
points(value~TminAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),], col=paste("#990000",a,sep=""), pch=19, cex=.5)
H <- lm(value~TminAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),])
M <- lm(value~TminAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),])
L <- lm(value~TminAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),])
abline(H, col="#000099", lwd=3)
abline(M, col="#990099", lwd=3)
abline(L, col="#CC0000", lwd=3)
mtext(text="c)", 3, adj=0.05, font=2, padj=1.5)




#----------------------------------------------------------------------
#========= Now with one climate variable but elevations in panels =====
#                          older Figure, opacity needs tweeked 
#       ALSO, DOES NOT HOLD YLIM CONSTANT!!!
quartz(width=7, height=5)
par(mfrow=c(1,3), oma=c(0,5,1,2), mar=c(7,0,3,0))
#xyplot(value~TminAn + TaveAn + PPTGS | Elev, groups=Comp, data=rw.clim.allraw, pch=19, type =c("p","r"),main=paste("PIPO-CMDGS"),lwd=3, drop.unused.levels=FALSE, col=c(rgb(.8,0,0,0.2), rgb(0,0,.6,0.3)))
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),], col="#99000030", pch=19, cex=.5, ylab="Ring Width Index", cex.axis=1.3, main="Low Elevation", xlab="")
H <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),])
M <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),])
L <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),])
#abline(H, col="gray", lwd=3 )
#abline(M, col="gray", lwd=3)
abline(L, col="#CC0000", lwd=3)

plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),], col="#99009955", pch=19, cex=.5, xlab="", cex.axis=1.3, font.lab=2, main= "Mid Elevation", yaxt="n")
abline(M, col="#990099", lwd=3)
#abline(H, col="gray", lwd=3, lty=2) # "#00009988"
#abline(L, col="gray", lwd=3, lty=2) # "#CC000088"

plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], col="#00009955", pch=19, cex=.5, main="High Elevation", xlab="", yaxt="n", cex.axis=1.3)
#abline(L, col="gray", lwd=3) #"#CC000088"
#abline(M, col="gray", lwd=3) # ""#99009999"
abline(H, col="#000099", lwd=3)
mtext(text="Ring Width Index", 2, outer=T, padj=-3, font=2, adj=.6)
mtext(text="Climate Moisture Deficit (PET-P)", 1, outer=T, padj=-3, font=2)



#######################################################################
#---------------- CMD in single panel for UW internal grant -----------
#                        One small plot, with just CMD
#----------------------------------------------------------------------
quartz(width=2.6, height=3.3)
a <-"06"
par(mar=c(3,3,1,1))
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], col=paste("#000099",a,sep=""),
     xlim=c(200,700),pch=19, cex=.5,
     ylab="Ring Width Index", xlab="Climate Moisture Deficit", cex.lab=1,
     font.lab=2, cex.axis=1, mgp=c(1.5,.4,0), tcl=-0.3)
points(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),], col=paste("#990099",a,sep=""), pch=19, cex=.5)
points(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),], col=paste("#990000",a,sep=""), pch=19, cex=.5)
H <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),])
M <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),])
L <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),])
abline(H, col="#000099", lwd=3)
abline(M, col="#990099", lwd=3)
abline(L, col="#CC0000", lwd=3)

#par(font=2) #can't make bold in legend function
legend("topright", legend=c("High", "Mid", "Low"), lwd=3, col=c("#000099","#990099","#CC0000"),bty="n", ncol=1, cex=0.9)
#par(font=1)



####################################################################################
########## Plotting to visualize the necessity of individual random effects ##############
#---------  for visualizing slopes of each individual from random effects -----------
#                  May not be working for lower elevations, and needs alpha work
#     -but take home point is it's reasonable for individual random effects to be a whash


testH <- lmer(value~CMDGS + (CMDGS|Tree), data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),])
testM <- lmer(value~CMDGS + (CMDGS|Tree), data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),])
testL <- lmer(value~CMDGS + (CMDGS|Tree), data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),])


#----------- Plotting individual lines underneath tree mean on plot from above -----------
quartz(width=7, height=5)
par(mfrow=c(1,3), oma=c(0,5,1,2), mar=c(7,0,3,0))
# Low ELEVATION
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),], col="#99000030", pch=19, cex=.5, ylab="Ring Width Index", cex.axis=1.3, main="Low Elevation", xlab="")
for(n in 1:nrow(coef(testL)$Tree)){
  abline(a = coef(testL)$Tree[n,1], b= coef(testL)$Tree[n,2], lwd=2, col="gray")
  print(n)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testL))[,1], lwd=3, col="#CC0000")

#MID ELEVATION
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),], col="#99009955", pch=19, cex=.5, xlab="", cex.axis=1.3, font.lab=2, main= "Mid Elevation", yaxt="n")
for(n in 1:nrow(coef(testM)$Tree)){
  abline(a = coef(testM)$Tree[n,1], b= coef(testM)$Tree[n,2], lwd=2, col="gray")
  print(n)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testM))[,1], lwd=3, col="#990099")

# HIGH ELEVATION
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], col="#00009955", pch=19, cex=.5, main="High Elevation", xlab="", yaxt="n", cex.axis=1.3)
for(n in 1:nrow(coef(testH)$Tree)){
  abline(a = coef(testH)$Tree[n,1], b= coef(testH)$Tree[n,2], lwd=2, col="gray")
  print(n)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testH))[,1], lwd=3, col="#000099")
mtext(text="Ring Width Index", 2, outer=T, padj=-3, font=2, adj=.6)
mtext(text="Climate Moisture Deficit (PET-P)", 1, outer=T, padj=-3, font=2)

### Seems to illustrate pretty well that there's almost no individual variation with this random effects structure

#------------------------------------------------
#         Plotting the same as above, but with individual regression lines
#         Rather than the shrinkage induced random effects lines
#--------------------
quartz(width=7, height=5)
par(mfrow=c(1,3), oma=c(0,5,1,2), mar=c(7,0,3,0))
# Low ELEVATION
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),], col="#99000025", pch=19, cex=.5, ylab="Ring Width Index", cex.axis=1.3, main="Low Elevation", xlab="")
for(n in 1:nrow(coef(testL)$Tree)){
  tree <- rownames(coef(testL)$Tree)[n]
  tmp <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L" & rw.clim.allraw$Tree==tree),])
  abline(tmp, lwd=1, col="#888888")
  print(tree)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testL))[,1], lwd=3, col="#CC0000")

#MID ELEVATION
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),], col="#99009925", pch=19, cex=.5, xlab="", cex.axis=1.3, font.lab=2, main= "Mid Elevation", yaxt="n")

for(n in 1:nrow(coef(testM)$Tree)){
  tree <- rownames(coef(testM)$Tree)[n]
  tmp <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M" & rw.clim.allraw$Tree==tree),])
  abline(tmp, lwd=1, col="#888888")
  print(tree)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testM))[,1], lwd=3, col="#990099")

# HIGH ELEVATION
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], col="#00009925", pch=19, cex=.5, main="High Elevation", xlab="", yaxt="n", cex.axis=1.3)
for(n in 1:nrow(coef(testH)$Tree)){
  tree <- rownames(coef(testH)$Tree)[n]
  tmp <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H" & rw.clim.allraw$Tree==tree),])
  abline(tmp, lwd=1, col="#888888")
  print(tree)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testH))[,1], lwd=3, col="#000099")
mtext(text="Ring Width Index", 2, outer=T, padj=-3, font=2, adj=.6)
mtext(text="Climate Moisture Deficit (PET-P)", 1, outer=T, padj=-3, font=2)

### Seems to illustrate pretty well that there's not much individual variation
# PARTICULARLY AT LOW ELEVATION
# now make a small plot with the slopes of the elevation regressions to show in corner
quartz(width=2, height=2)
par(mar=c(.5,2.5,.5,.2))
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], type="n", xaxt="n",yaxt="n", xlab="", ylab="", cex.axis=.8,bty="n", ylim=c(0.5,1.5))#, col="#00009925", pch=19, cex=.5, main="High Elevation", xlab="", yaxt="n", cex.axis=1.3)
axis(2,at=c(0.5,1,1.5),labels=c("0.5", "1.0", "1.5"))
abline(coef=coef(summary(testH))[,1], lwd=5, col="#0000DD")
abline(coef=coef(summary(testM))[,1], lwd=5, col="#DD00DD")
abline(coef=coef(summary(testL))[,1], lwd=5, col="#CC0000")





#=====================================================================================
########## Making six panel plot, elevations above, climate variables below ###########
#----------------------final plot for R Graphics class      -----------------------
#=====================================================================================

# makes a very usable plot. I edited it in powerpoint and made .pdf and .png formats
# from powerpoint. however, be careful moving powerpoint slides to Windows computers
# Should save as something other than a pdf and import again if I want a powerpoint slide

quartz(width=7, height=6.5)
par(mfrow=c(2,3), oma=c(0,5,2.5,2), mar=c(5,0,0,0))
textcex <- .8
a <-"09"
indivs <- "#B9B9B9"
ylims <- c(0,2.2)
# Low ELEVATION
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),], ylim=ylims, col=paste("#990000",a,sep=""), pch=19, cex=.5, ylab="Ring Width Index", cex.axis=1.3, xlab="")
mtext(text="Low Elevation",side=3,line=-1.5, cex=textcex, font=1)
for(n in 1:nrow(coef(testL)$Tree)){
  tree <- rownames(coef(testL)$Tree)[n]
  tmp <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L" & rw.clim.allraw$Tree==tree),])
  abline(tmp, lwd=1, col=indivs)
  print(tree)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testL))[,1], lwd=3, col="#CC0000")
mtext(text="a)", 3, adj=0.05, font=2, padj=1.5)

#MID ELEVATION
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),], ylim=ylims, col=paste("#990099",a,sep=""), pch=19, cex=.5, xlab="", cex.axis=1.3, font.lab=2, yaxt="n")
mtext(text="Mid Elevation",side=3,line=-1.5, cex=textcex, font=1)
for(n in 1:nrow(coef(testM)$Tree)){
  tree <- rownames(coef(testM)$Tree)[n]
  tmp <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M" & rw.clim.allraw$Tree==tree),])
  abline(tmp, lwd=1, col=indivs)
  print(tree)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testM))[,1], lwd=3, col="#990099")
mtext(text="Climate Moisture Deficit (PET-P, mm water deficit)", 3, outer=T, padj=0, font=2, cex=1)
mtext(text="b)", 3, adj=0.05, font=2, padj=1.5)

# HIGH ELEVATION
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], ylim=ylims, col=paste("#000099",a,sep=""), pch=19, cex=.5, xlab="", yaxt="n", cex.axis=1.3)
mtext(text="High Elevation",side=3,line=-1.5, cex=textcex, font=1)
for(n in 1:nrow(coef(testH)$Tree)){
  tree <- rownames(coef(testH)$Tree)[n]
  tmp <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H" & rw.clim.allraw$Tree==tree),])
  abline(tmp, lwd=1, col=indivs)
  print(tree)
}
mtext(text="c)", 3, adj=0.05, font=2, padj=1.5)
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testH))[,1], lwd=3, col="#000099")
mtext(text="Ring Width Index", 2, outer=T, padj=-3, font=2, adj=.85)


#------- now lower row with different climate variables ------

a<-"09" # choose an alpha value for opacity
#xyplot(value~TminAn + TaveAn + PPTGS | Elev, groups=Comp, data=rw.clim.allraw, pch=19, type =c("p","r"),main=paste("PIPO-CMDGS"),lwd=3, drop.unused.levels=FALSE, col=c(rgb(.8,0,0,0.2), rgb(0,0,.6,0.3)))
plot(value~TmaxAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], ylim=ylims, col=paste("#000099",a,sep=""), xlim=c(22,28.3), pch=19, cex=.5, ylab="Ring Width Index", xlab="Max Temp (°C)", cex.lab=1.4, font.lab=2, cex.axis=1.3)
points(value~TmaxAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),], col=paste("#990099",a,sep=""), pch=19, cex=.5)
points(value~TmaxAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),], col=paste("#990000",a,sep=""), pch=19, cex=.5)
H <- lm(value~TmaxAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),])
M <- lm(value~TmaxAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),])
L <- lm(value~TmaxAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),])
abline(H, col="#000099", lwd=3)
abline(M, col="#990099", lwd=3)
abline(L, col="#CC0000", lwd=3)
mtext(text="Ring Width Index", 2, outer=T, padj=-3, font=2, adj=.25)
mtext(text="d)", 3, adj=0.05, font=2, padj=1.5)
legend("topright", legend=c("High", "Mid", "Low"), lwd=3, col=c("#000099","#990099","#CC0000"),bty="n", cex=1)

plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], ylim=ylims, col=paste("#000099",a,sep=""), xlim=c(220,680),pch=19, cex=.5, yaxt="n", xlab="CMD (mm deficit)", cex.lab=1.4, font.lab=2, cex.axis=1.3)
points(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),], col=paste("#990099",a,sep=""), pch=19, cex=.5)
points(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),], col=paste("#990000",a,sep=""), pch=19, cex=.5)
H <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),])
M <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),])
L <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),])
abline(H, col="#000099", lwd=3)
abline(M, col="#990099", lwd=3)
abline(L, col="#CC0000", lwd=3)
mtext(text="e)", 3, adj=0.05, font=2, padj=1.5)

#par(font=2) #can't make bold in legend function
#legend(x=450, y=2.33, xpd=NA, xjust=0.5, legend=c("High Elev", "Mid Elev", "Low Elev"), lwd=3, col=c("#000099","#990099","#CC0000"),bty="n", ncol=3, cex=1.2)
#legend("topright", xjust=0.5, legend=c("High Elev", "Mid Elev", "Low Elev"), lwd=3, col=c("#000099","#990099","#CC0000"),bty="n", ncol=0, cex=1)
#par(font=1)
# x=280, y=2.0 to get it in the plot, the xjust gives me center justification on x value, default seems to be left

plot(value~TminAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], ylim=ylims, col=paste("#000099",a,sep=""), pch=19, cex=.5, yaxt="n", xlab="Min Temp (°C)", cex.lab=1.4, font.lab=2, cex.axis=1.3)
points(value~TminAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),], col=paste("#990099",a,sep=""), pch=19, cex=.5)
points(value~TminAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),], col=paste("#990000",a,sep=""), pch=19, cex=.5)
H <- lm(value~TminAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),])
M <- lm(value~TminAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),])
L <- lm(value~TminAn, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),])
abline(H, col="#000099", lwd=3)
abline(M, col="#990099", lwd=3)
abline(L, col="#CC0000", lwd=3)
mtext(text="f)", 3, adj=0.05, font=2, padj=1.5)




#################### Working with ASPENS ################

#### NOTE: THINGS DON'T WORK SUPER WELL FROM HERE ON DOWN

rw.clim.allelev <- read.csv(file="CO-POTR-allRWIClim.csv")
rw.clim.allraw <- read.csv(file="CO-POTR-allRWI_rawClim.csv")

####################################################################################
########## Plotting to visualize the necessity of individual random effects ##############
#---------  for visualizing slopes of each individual from random effects -----------
#                  May not be working for lower elevations, and needs alpha work
#     -but take home point is it's reasonable for individual random effects to be a whash


testH <- lmer(value~CMDGS + (CMDGS|Tree), data=rw.clim.allelev[which(rw.clim.allelev$Elev=="H"),])
testM <- lmer(value~CMDGS + (CMDGS|Tree), data=rw.clim.allelev[which(rw.clim.allelev$Elev=="M"),])
testL <- lmer(value~CMDGS + (CMDGS|Tree), data=rw.clim.allelev[which(rw.clim.allelev$Elev=="L"),])


#----------- Plotting individual lines underneath tree mean on plot from above -----------
quartz(width=7, height=5)
par(mfrow=c(1,3), oma=c(0,5,1,2), mar=c(7,0,3,0))
# Low ELEVATION
plot(value~CMDGS, data=rw.clim.allelev[which(rw.clim.allelev$Elev=="L"),], col="#99000030", pch=19, cex=.5, ylab="Ring Width Index", cex.axis=1.3, main="Low Elevation", xlab="")
for(n in 1:nrow(coef(testL)$Tree)){
  abline(a = coef(testL)$Tree[n,1], b= coef(testL)$Tree[n,2], lwd=2, col="gray")
  print(n)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testL))[,1], lwd=3, col="#CC0000")

#MID ELEVATION
plot(value~CMDGS, data=rw.clim.allelev[which(rw.clim.allelev$Elev=="M"),], col="#99009955", pch=19, cex=.5, xlab="", cex.axis=1.3, font.lab=2, main= "Mid Elevation", yaxt="n")
for(n in 1:nrow(coef(testM)$Tree)){
  abline(a = coef(testM)$Tree[n,1], b= coef(testM)$Tree[n,2], lwd=2, col="gray")
  print(n)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testM))[,1], lwd=3, col="#990099")

# HIGH ELEVATION
plot(value~CMDGS, data=rw.clim.allelev[which(rw.clim.allelev$Elev=="H"),], col="#00009955", pch=19, cex=.5, main="High Elevation", xlab="", yaxt="n", cex.axis=1.3)
for(n in 1:nrow(coef(testH)$Tree)){
  abline(a = coef(testH)$Tree[n,1], b= coef(testH)$Tree[n,2], lwd=2, col="gray")
  print(n)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testH))[,1], lwd=3, col="#000099")
mtext(text="Ring Width Index", 2, outer=T, padj=-3, font=2, adj=.6)
mtext(text="Climate Moisture Deficit (PET-P)", 1, outer=T, padj=-3, font=2)


#------------------------------------------------
#         Plotting the same as above, but with individual regression lines
#         Rather than the shrinkage induced random effects lines
#--------------------
quartz(width=7, height=5)
par(mfrow=c(1,3), oma=c(0,5,1,2), mar=c(7,0,3,0))
# Low ELEVATION
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L"),], col="#99000025", pch=19, cex=.5, ylab="Ring Width Index", cex.axis=1.3, main="Low Elevation", xlab="")
for(n in 1:nrow(coef(testL)$Tree)){
  tree <- rownames(coef(testL)$Tree)[n]
  tmp <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="L" & rw.clim.allraw$Tree==tree),])
  abline(tmp, lwd=1, col="#888888")
  print(tree)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testL))[,1], lwd=3, col="#CC0000")

#MID ELEVATION
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M"),], col="#99009925", pch=19, cex=.5, xlab="", cex.axis=1.3, font.lab=2, main= "Mid Elevation", yaxt="n")

for(n in 1:nrow(coef(testM)$Tree)){
  tree <- rownames(coef(testM)$Tree)[n]
  tmp <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="M" & rw.clim.allraw$Tree==tree),])
  abline(tmp, lwd=1, col="#888888")
  print(tree)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testM))[,1], lwd=3, col="#990099")

# HIGH ELEVATION
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], col="#00009925", pch=19, cex=.5, main="High Elevation", xlab="", yaxt="n", cex.axis=1.3)
for(n in 1:nrow(coef(testH)$Tree)){
  tree <- rownames(coef(testH)$Tree)[n]
  tmp <- lm(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H" & rw.clim.allraw$Tree==tree),])
  abline(tmp, lwd=1, col="#888888")
  print(tree)
}
#abline(H, col="#CC0000", lwd=3) # this is almost identicial to testH, because random effects don't seem to matter
abline(coef=coef(summary(testH))[,1], lwd=3, col="#000099")
mtext(text="Ring Width Index", 2, outer=T, padj=-3, font=2, adj=.6)
mtext(text="Climate Moisture Deficit (PET-P)", 1, outer=T, padj=-3, font=2)

### Seems to illustrate pretty well that there's not much individual variation
# PARTICULARLY AT LOW ELEVATION
# now make a small plot with the slopes of the elevation regressions to show in corner
quartz(width=2, height=2)
par(mar=c(.5,2.5,.5,.2))
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], type="n", xaxt="n",yaxt="n", xlab="", ylab="", cex.axis=.8,bty="n", ylim=c(0.5,1.5))#, col="#00009925", pch=19, cex=.5, main="High Elevation", xlab="", yaxt="n", cex.axis=1.3)
axis(2,at=c(0.5,1,1.5),labels=c("0.5", "1.0", "1.5"))
abline(coef=coef(summary(testH))[,1], lwd=5, col="#0000DD")
abline(coef=coef(summary(testM))[,1], lwd=5, col="#DD00DD")
abline(coef=coef(summary(testL))[,1], lwd=5, col="#CC0000")






