######################################
#       Two panel plot of PIPO
#     CMI vs RWI and BAI10yr vs elevxcomp
#             for: DDIG proposal 2014
#######################################
# Requires:
#       "PIPO-BAI10yr.csv" and  
#       "CO-PIPO-allRWI_rawClim.csv" 
# and makes a one panel dot w/ elevation regression plot (panel A)
# and another panel for elevationxcomp BAI10yr plot in nice figure

########### Code lifted from Fig_PIPO_RW_v_Clim.r

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


############# Import BAI data #######
tree_info <- read.csv("PIPO-BAI10yr.csv")
levels(tree_info$Band) <- list(L="L", M = "M", H = "H")
levels(tree_info$Comp) <- list(H = "H", L="L")


#######################################################################
#---------------- CMD in single panel for UW internal grant -----------
#                        One small plot, with just CMD
#----------------------------------------------------------------------
quartz(width=5.5, height=3.3)
a <-"06"
par(mfrow=c(1,2), mar=c(3,3,1.1,1), mgp=c(1.5,.4,0), cex.lab=1.05, cex.axis=.8, font.lab=1, tcl=-0.3)
plot(value~CMDGS, data=rw.clim.allraw[which(rw.clim.allraw$Elev=="H"),], col=paste("#000099",a,sep=""),
     xlim=c(200,700),pch=19, cex=.5,
     ylab="Ring Width Index", xlab="Climate Moisture Deficit")
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



test <- expression("10yr Mean Basal Area Increment (m" * m^2 * ")")

######### code lifted from Code-Competition_analisis11_11_13

boxplot(BAI10yr~Comp*Band, data=tree_info, col=c("#CC00005E", "#CC0000EE","#99009955","#990099DD","#00009940","#000099CC" ),boxwex=.7, outwex=0
        , at = c(1,2,4,5,7,8)
        , ylab = test
        , xlab = "Elevation Band"
        , xaxt = "n"
        , main = ""
        , mgp=c(1.5,.4,0)
)
axis(1,at=c(1.5,4.5,7.5), labels=c("Low", "Mid", "High"))
legend("topleft",legend=c("Competitive", "Noncompetitive"),fill=c("#00009940","#000099CC"), bty="n", cex=.9)
