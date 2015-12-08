###################################
## Combining summary stats from Synchrony and BAI-Elev-Competition ###
###################################


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# this code takes the outputs of the Synchrony and BAI/ElevComp analysis (plus hopefully future analyses)
# and mashes them together to create a master dataframe for plotting/synthesis

# **this code is not stand-alone**

# last updated: 12/7/15

# Note: for the moment, you must run:
# Cd-Synchrony.R
# Cd-BAI-vs-ElevComp.R
# before running this code in order to store BAI models and Synchrony sumstats dfs
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




####### Best BAI and BAIcv estimates from BAI-vs-ElevComp code: ######
# this code extracts the model coefficients from the best fit BAI model determined in Cd-BAI0-vs-ElevComp.R
# and stores mean BAI (+se) and BAI CV (+se) estimates from these models in the sumstats dataframe for each species

### creating xvales that work with a standard barplot so that bars and points plot in the appropriate location
xvals <- c(0.6666,1.66666,2.6666)

### CO_PIPO ##
meanmod <- CO_PIPOcompmod # get the best fit BAI ~ Comp + Elev model for this species
CVmod <- CO_PIPOcompCV # get the best BAIcv ~ Comp + Elev model for this species
dataz <- CO_PIPOinfo # get the underlying data for this species so I can update the model
meanmod1 <- update(meanmod, .~.-1) # update the model removing intercept to make coefficients easy to interpret
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2]) # extract coefficients and se's
meangr$xval <- xvals # add the xvals purely for plotting purposes
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
CO_PIPOsumstats1 <- cbind(CO_PIPOsumstats, meangr, CVgr)
colnames(CO_PIPOsumstats1)[9:13] <- c("BAI.mod.mean", "BAI.mod.se", "xval", "BAIcv.mod.mean", "BAIcv.mod.se")

### CO_POTR ##
meanmod <- CO_POTRcompmod
CVmod <- CO_POTRcompCV
dataz <- CO_POTRinfo
meanmod1 <- update(CO_POTRcompmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
CO_POTRsumstats1 <- cbind(CO_POTRsumstats, meangr, CVgr)
colnames(CO_POTRsumstats1)[9:13] <- c("BAI.mod.mean", "BAI.mod.se", "xval", "BAIcv.mod.mean", "BAIcv.mod.se")


### CO_ABLA ##
meanmod <- CO_ABLAcompmod
CVmod <- CO_ABLAcompCV
dataz <- CO_ABLAinfo
meanmod1 <- update(meanmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
CO_ABLAsumstats1 <- cbind(CO_ABLAsumstats, meangr, CVgr)
colnames(CO_ABLAsumstats1)[9:13] <- c("BAI.mod.mean", "BAI.mod.se", "xval", "BAIcv.mod.mean", "BAIcv.mod.se")



#### MONTANA _________________________________
### MT-TSHE ##
meanmod <- MT_TSHEcompmod
CVmod <- MT_TSHEcompCV
dataz <- MT_TSHEinfo
meanmod1 <- update(meanmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
MT_TSHEsumstats1 <- cbind(MT_TSHEsumstats, meangr, CVgr)
colnames(MT_TSHEsumstats1)[9:13] <- c("BAI.mod.mean", "BAI.mod.se", "xval", "BAIcv.mod.mean", "BAIcv.mod.se")


### MT-PSME ##
meanmod <- MT_PSMEcompmod
CVmod <- MT_PSMEcompCV
dataz <- MT_PSMEinfo
meanmod1 <- update(CO_POTRcompmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
MT_PSMEsumstats1 <- cbind(MT_PSMEsumstats, meangr, CVgr)
colnames(MT_PSMEsumstats1)[9:13] <- c("BAI.mod.mean", "BAI.mod.se", "xval", "BAIcv.mod.mean", "BAIcv.mod.se")


### MT_ABLA ##
meanmod <- MT_ABLAcompmod
CVmod <- MT_ABLAcompCV
dataz <- MT_ABLAinfo
meanmod1 <- update(meanmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
MT_ABLAsumstats1 <- cbind(MT_ABLAsumstats, meangr, CVgr)
colnames(MT_ABLAsumstats1)[9:13] <- c("BAI.mod.mean", "BAI.mod.se", "xval", "BAIcv.mod.mean", "BAIcv.mod.se")



#### WASHINGTON _________________________________
### WA-TSHE ##
meanmod <- WA_TSHEcompmod
CVmod <- WA_TSHEcompCV
dataz <- WA_TSHEinfo
meanmod1 <- update(meanmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:2,1:2])
meangr$xval <- c(1.66666,2.66666)
CVgr <- data.frame(summary(CVmod1)$coefficients[1:2,1:2])
WA_TSHEsumstats1 <- cbind(WA_TSHEsumstats, meangr, CVgr)
colnames(WA_TSHEsumstats1)[9:13] <- c("BAI.mod.mean", "BAI.mod.se", "xval", "BAIcv.mod.mean", "BAIcv.mod.se")


### MT-PSME ##
meanmod <- WA_PSMEcompmod
CVmod <- WA_PSMEcompCV
dataz <- WA_PSMEinfo
meanmod1 <- update(CO_POTRcompmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
WA_PSMEsumstats1 <- cbind(WA_PSMEsumstats, meangr, CVgr)
colnames(WA_PSMEsumstats1)[9:13] <- c("BAI.mod.mean", "BAI.mod.se", "xval", "BAIcv.mod.mean", "BAIcv.mod.se")

### WA_ABLA ##
meanmod <- WA_ABLAcompmod
CVmod <- WA_ABLAcompCV
dataz <- WA_ABLAinfo
meanmod1 <- update(meanmod, .~.-1)
CVmod1 <- update(CVmod, .~.-1)
meangr <- data.frame(summary(meanmod1)$coefficients[1:3,1:2])
meangr$xval <- xvals
CVgr <- data.frame(summary(CVmod1)$coefficients[1:3,1:2])
WA_ABLAsumstats1 <- cbind(WA_ABLAsumstats, meangr, CVgr)
colnames(WA_ABLAsumstats1)[9:13] <- c("BAI.mod.mean", "BAI.mod.se", "xval", "BAIcv.mod.mean", "BAIcv.mod.se")


#______ Combine all species sumstats into one large dataframe ___________
Sumstatsall <- rbind(CO_PIPOsumstats1,CO_POTRsumstats1, CO_ABLAsumstats1,
                     MT_TSHEsumstats1,MT_PSMEsumstats1, MT_ABLAsumstats1,
                     WA_TSHEsumstats1,WA_PSMEsumstats1, WA_ABLAsumstats1)

Sumstatsall$Species <- as.factor(c(rep("PIPO", times=3), rep("POTR", times=3), rep("ABLA",times=3),
                         rep("TSHE", times=3), rep("PSME", times=3), rep("ABLA", times=3),
                         rep("TSHE", times=2), rep("PSME", times=3), rep("ABLA", times=3)))
Sumstatsall$State <- as.factor(c(rep("CO", times=9), rep("MT", times=9), rep("WA", times=8)))


######### END: BAI and BAIcv estimates ###################################





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##### ******Plotting various synthesis plots of BAI, BAIcv, synchrony, etc.********** ##############


###### Plotting BAI against Synchrony, all in one plot ######
pal <- brewer.pal(n = 7, name="Dark2")
palette(pal)
plot(corr.mod.mean~BAI.mod.mean, Sumstatsall, pch=15+as.numeric(State), col=as.factor(Species), cex=1.4, typ="n", ylab="growth synchrony", xlab="mean BAI (cm^3)")
for(i in unique(Sumstatsall$State)){
  for (j in unique(Sumstatsall$Species)){
    points(corr.mod.mean~BAI.mod.mean, Sumstatsall, subset=Species==j & State==i, pch=15+as.numeric(State), col=as.factor(Species), cex=1.4, type="b")
  }
}


####### Plotting BAI agains BAI CV all in one plot ########
plot(BAIcv.mod.mean~BAI.mod.mean, Sumstatsall, pch=15+as.numeric(State), col=as.factor(Species), cex=1.4, typ="n", ylab="BAI CV", xlab="mean BAI (cm^3)")
for(i in unique(Sumstatsall$State)){
  for (j in unique(Sumstatsall$Species)){
    points(BAIcv.mod.mean~BAI.mod.mean, Sumstatsall, subset=Species==j & State==i, pch=15+as.numeric(State), col=as.factor(Species), cex=1.4, type="b")
  }
}


####### Plotting BAIcv agains Synchrony all in one plot ######
plot(BAIcv.mod.mean~corr.mod.mean, Sumstatsall, pch=15+as.numeric(State), col=as.factor(Species), cex=1.4, typ="n", ylab="BAI CV", xlab="mean BAI (cm^3)")
for(i in unique(Sumstatsall$State)){
  for (j in unique(Sumstatsall$Species)){
    points(BAIcv.mod.mean~corr.mod.mean, Sumstatsall, subset=Species==j & State==i, pch=15+as.numeric(State), col=as.factor(Species), cex=1.4, type="b")
  }
}


###### Plotting BAI~Synchrony and BAI~BAICV for each state ####
quartz(width=6, height=4.3)
par(mfrow=c(2,3), mar=c(0,3,0,0), oma=c(4,2.2,2,1))

pal <- brewer.pal(n = 7, name="Dark2")
palette(pal)

for(i in unique(Sumstatsall$State)){
  plot(corr.mod.mean~BAI.mod.mean, Sumstatsall,subset=State==i, pch=15+as.numeric(State), col=as.factor(Species), cex=1.4, typ="n", ylab="growth synchrony", xaxt='n')
  for (j in unique(Sumstatsall$Species)){
    points(corr.mod.mean~BAI.mod.mean, Sumstatsall, subset=Species==j & State==i, pch=15+as.numeric(State), col=as.factor(Species), cex=1.4, type="b")
  }
}


####### Plotting BAI agains BAI CV all in one plot ########
for(i in unique(Sumstatsall$State)){
  plot(BAIcv.mod.mean~BAI.mod.mean, Sumstatsall, subset=State==i, pch=15+as.numeric(State), col=as.factor(Species), cex=1.4, typ="n", ylab="BAI CV", xlab="mean BAI (cm^3)")
  for (j in unique(Sumstatsall$Species)){
    points(BAIcv.mod.mean~BAI.mod.mean, Sumstatsall, subset=Species==j & State==i, pch=15+as.numeric(State), col=as.factor(Species), cex=1.4, type="b")
  }
}

mtext(text = "mean BAI (cm^3)", side=1, outer = T, line=2.5)
mtext(text = "Growth Synchrony", side=2, outer=T, adj=0.85)
mtext(text = "BAI CV", side=2, outer=T, adj=0.15)
mtext(text = "Colorado", side=3, outer=T, adj=0.15, cex=1)
mtext(text = "Montana", side=3, outer=T, adj=0.55, cex=1)
mtext(text = "Washington", side=3, outer=T, adj=0.95, cex=1)
legend("bottomright", bty="n", legend = c('ABLA',"PIPO","POTR","PSME","TSHE")
       , col=pal[1:5], lwd=1, pch=16, cex = .8)



############ 3 panel figures of Synchrony ##########
quartz(width=6, height=3)
par(mfrow=c(1,3), mar=c(1,3,1,0), oma=c(0,2,0,1))
xvals <- c(0.6666,1.66666,2.6666)
colors <- brewer.pal(n=5, "Set1")[c(4,5,3)] 

varname <- "Meancorr"
errorname <- "SDcorr"

plot(get(varname)~xval, data=Sumstatsall,subset=Species=="PIPO" & State=="CO", xlim=c(0,3.333), bty="n", xaxt="n", ylab="Growth Synchrony", col=colors[1], ylim=c(0,1), type="b" ) 
error_bars(xvar ="xval", yvar=varname,upper=errorname,errordata = Sumstatsall[which(Sumstatsall$State=="CO"&Sumstatsall$Species=="PIPO"),],length = 0, lwd=2, col=colors[1]) 

plot(get(varname)~xval, data=Sumstatsall,subset=Species=="POTR" & State=="CO", xlim=c(0,3.333), bty="n", xaxt="n", ylab="Growth Synchrony", col=colors[2], ylim=c(0,1), type="b" ) 
error_bars(xvar ="xval", yvar=varname,upper=errorname,errordata = Sumstatsall[which(Sumstatsall$State=="CO"&Sumstatsall$Species=="POTR"),],length = 0, lwd=2, col=colors[2]) 

plot(get(varname)~xval, data=Sumstatsall,subset=Species=="ABLA" & State=="CO", xlim=c(0,3.333), bty="n", xaxt="n", ylab="Growth Synchrony", col=colors[3], ylim=c(0,1), type="b" ) 
error_bars(xvar ="xval", yvar=varname,upper=errorname,errordata = Sumstatsall[which(Sumstatsall$State=="CO"&Sumstatsall$Species=="ABLA"),],length = 0, lwd=2, col=colors[3]) 



#### Montana 
quartz(width=6, height=3)
par(mfrow=c(1,3), mar=c(1,3,1,0), oma=c(0,2,0,1))
xvals <- c(0.6666,1.66666,2.6666)
colors <- brewer.pal(n=3, "Set1") 

varname <- "Meancorr"
errorname <- "SDcorr"

plot(get(varname)~xval, data=Sumstatsall,subset=Species=="TSHE" & State=="MT", xlim=c(0,3.333), bty="n", xaxt="n", ylab="Growth Synchrony", col=colors[1], ylim=c(0,1), type="b" ) 
error_bars(xvar ="xval", yvar=varname,upper=errorname,errordata = Sumstatsall[which(Sumstatsall$State=="MT"&Sumstatsall$Species=="TSHE"),],length = 0, lwd=2, col=colors[1]) 

plot(get(varname)~xval, data=Sumstatsall,subset=Species=="PSME" & State=="MT", xlim=c(0,3.333), bty="n", xaxt="n", ylab="Growth Synchrony", col=colors[2], ylim=c(0,1), type="b" ) 
error_bars(xvar ="xval", yvar=varname,upper=errorname,errordata = Sumstatsall[which(Sumstatsall$State=="MT"&Sumstatsall$Species=="PSME"),],length = 0, lwd=2, col=colors[2]) 

plot(get(varname)~xval, data=Sumstatsall,subset=Species=="ABLA" & State=="MT", xlim=c(0,3.333), bty="n", xaxt="n", ylab="Growth Synchrony", col=colors[3], ylim=c(0,1), type="b" ) 
error_bars(xvar ="xval", yvar=varname,upper=errorname,errordata = Sumstatsall[which(Sumstatsall$State=="MT"&Sumstatsall$Species=="ABLA"),],length = 0, lwd=2, col=colors[3]) 



#### Washington
quartz(width=6, height=3)
par(mfrow=c(1,3), mar=c(1,3,1,0), oma=c(0,2,0,1))
xvals <- c(0.6666,1.66666,2.6666)
colors <- brewer.pal(n=3, "Set1") 

varname <- "Meancorr"
errorname <- "SDcorr"

plot(get(varname)~xval, data=Sumstatsall,subset=Species=="TSHE" & State=="WA", xlim=c(0,3.333), bty="n", xaxt="n", ylab="Growth Synchrony", col=colors[1], ylim=c(0,1), type="b" ) 
error_bars(xvar ="xval", yvar=varname,upper=errorname,errordata = Sumstatsall[which(Sumstatsall$State=="WA"&Sumstatsall$Species=="TSHE"),],length = 0, lwd=2, col=colors[1]) 

plot(get(varname)~xval, data=Sumstatsall,subset=Species=="PSME" & State=="WA", xlim=c(0,3.333), bty="n", xaxt="n", ylab="Growth Synchrony", col=colors[2], ylim=c(0,1), type="b" ) 
error_bars(xvar ="xval", yvar=varname,upper=errorname,errordata = Sumstatsall[which(Sumstatsall$State=="WA"&Sumstatsall$Species=="PSME"),],length = 0, lwd=2, col=colors[2]) 

plot(get(varname)~xval, data=Sumstatsall,subset=Species=="ABLA" & State=="WA", xlim=c(0,3.333), bty="n", xaxt="n", ylab="Growth Synchrony", col=colors[3], ylim=c(0,1), type="b" ) 
error_bars(xvar ="xval", yvar=varname,upper=errorname,errordata = Sumstatsall[which(Sumstatsall$State=="WA"&Sumstatsall$Species=="ABLA"),],length = 0, lwd=2, col=colors[3]) 
