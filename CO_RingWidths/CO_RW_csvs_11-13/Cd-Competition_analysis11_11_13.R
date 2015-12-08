#=======================================================
#-------------------------------------------------------
#               Competition Analysis
#               Leander DL Anderegg
#                   11/11/13
#-------------------------------------------------------
#=======================================================
### Initially just written for PIPO, but expanded to POTR


rm(list=ls())
setwd("~/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13")

# set working directory to my scratch file in Dropbox

### load the core_average() function from  fn_core_average.R
source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/fn_core_average.R")
source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/fn_core_trunc.R")
source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/fn_sparkplotfn.R")

#===============================================================
############## PIPO Analysis ###################################
#===============================================================


#++++++++++++++++++++++++++++++++++++++++++++++
# NOTE: code between the ++++ markes resuses
  # nomenclature for different species. so must be rerun for each species
## low elevation
csv1 <- "CO-PIPO-L_RWdata11_7_13test2.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
treedata1 <- core_average(data=data1)
## mid elevation
csv2 <- "CO-PIPO-M_RWdata11_7_13.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "CO-PIPO-H_RWdata10_31_13.csv"
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3)

#------------- Getting 10yr Average Growth ---------------
tree10yr <- rbind(treedata1[,1:15], treedata2[1:15], treedata3[,1:15])
tenyrmean <- apply (X=tree10yr[,5:15],MARGIN=1, FUN=mean) # getting 10 yr average growth
Treat <- rep(c("C", "N"),times= length(tenyrmean)/2) 
    #can't think of a better way to get competitive treatment right now
    # but seems to work because H is always first. switching to C and N to avoid confusion

comp_PIPO <- cbind(tree10yr[1:4], Treat, tenyrmean)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

test <- aov (tenyrmean~Treat*Elev, data = comp_PIPO)


# --------- Preliminary Boxplot --------------------
quartz(width=4, height=4.5)
boxplot(tenyrmean~Treat*Elev, data=comp_PIPO, col=c("grey", "white")
        , ylab = "Average Ring Width (mm)"
        , xlab = "Elevation"
        , xaxt = "n"
        , main = "PIPO 10yr average growth")
axis(1,at=c(1.5,3.5,5.5), labels=c("Low", "Med", "High"))
legend("topleft",legend=c("Competitive", "Noncompetitive"),fill=c("grey", "white"))

# interaction plot
interaction.plot(x.factor=comp_PIPO$Elev, trace.factor=comp_PIPO$Treat, response=comp_PIPO$tenyrmean)


###############################################
#         BAI calculations                    #
###############################################
# spreadsheet to get DBH and bark depth out of
COtreesnew <- read.csv("/Users/leeanderegg/Dropbox/treecoresranges(leander) (1)/CO_RingWidths/CO_Checks/CO_focaltree_3_6.csv", header=T)

# make tags for common reference
ncols <- min(c(ncol(treedata1), ncol(treedata2), ncol(treedata3)))
rwdata_all <- rbind(treedata1[,1:ncols], treedata2[1:ncols], treedata3[1:ncols])
rwdata_all$Tree <- as.character(rwdata_all$Tree)
####### Need to change when switching nomenclature
rwdata_all$Tree[which(rwdata_all$Tree=="1L")] <- "01L" #having problems with 1 vs 10,11,12 later in script, so changing to 01
rwdata_all$Tree[which(rwdata_all$Tree=="1H")] <- "01H"
rwdata_all$Tag <- paste(rwdata_all$Species, rwdata_all$Elev, rwdata_all$Tree, sep = "-")
  #made unique tree tag
rwdata_alltag <- rwdata_all[,c(1:4,ncol(rwdata_all),5:(ncol(rwdata_all)-1))]
  #reordered columns to get tag in right place
rwdata_alltag <- rwdata_alltag[order(rwdata_alltag$Tag),]

# make a treeDBH data frame with just the tag, DBH and bark
tree_info <- COtreesnew[which(COtreesnew$Tag %in% rwdata_alltagPIPO$Tag),]
#tree_info <- COtreesnew[match(rwdata)]
#MtreeDBH <- COtreesnew[which(COtreesnew$Tag %in% treedata2$Tag), c(6,12,13)]
#LtreeDBH <- COtreesnew[which(COtreesnew$Tag %in% treedata1$Tag), c(6,12,13)]
tree_info <- tree_info[order(tree_info$Tag),] # get in in the same order as treedata
tree_infoPIPO <- tree_info
### MAKE SURE THEY:RE IN THE RIGHT ORDER !!!!!!!!!!!
rwdata_alltag$Tag==tree_info$Tag
rwdata_alltagPIPO <- rwdata_alltag
################ the function ######################
  # requires: 
      # ringwidths = ring width data frame from core_average(), with $Tag
      # DBHbark = tree DBH and bark info, with $Tag
      # DBHcol and barkcol = column numbers where to find DBH and bark
basal_area_inc <- function(ringwidths, DBHbark, DBHcol=12, barkcol=13){
  #make a Tag for ringwidths that works with my DBH data frame
  #get the initial radius ()
  radinit <- ((DBHbark[,DBHcol])/2)-DBHbark[,barkcol] # in cm
  DBHbark$radinit <- radinit *10 #in mm
  # order them so that they're in the same order
  DBHbark <- DBHbark[order(DBHbark$Tag),]
  ringwidths <- ringwidths[order(ringwidths$Tag),]
  # now calculate BAI for each year
  radius <- matrix(nrow=nrow(ringwidths), ncol=ncol(ringwidths)-5)
  for (i in 1:(ncol(ringwidths)-5)){ #account for the first five non rw columns
    if (i == 1)
      radius[,i] <- DBHbark$radinit
    else
      radius[,i] <- radius[,(i-1)]-ringwidths[,5+i]
    #print(names(ringwidths)[i+5])
    # Haha! only have one negative sized tree! boo yah!
  }
  print(paste("Min radius", min(radius, na.rm=T)))
  # ok, I've got the radii, now just need to convert to BAI
  BAI <- matrix(nrow=dim(radius)[1], ncol=dim(radius)[2])
  for(i in 1:(ncol(radius)-1)){ # this will return a matrix with 1 fewer column than radius
    BAI[,i] <- (radius[,i]^2 *pi) - (radius[,i+1]^2 *pi)
  }
  colnames(BAI) <- colnames(ringwidths)[-c(1:5)]
  rownames(BAI) <- DBHbark$Tag
  return(BAI)  
  
}
####################################################
require(reshape)
# rwdata_alltag <- read.csv("PIPO_avgRingWidth_2-19-15.csv", header=T)
# tree_info <- read.csv("PIPO_TreeInfo.csv", header=T)


BAI <- basal_area_inc(rwdata_alltag, tree_info)
min(BAI, na.rm=T)
# that seemed to do it. got a tree that's 8mm negative, but will just roll with it

#BAIn <- t(BAI)


#BAIplot <- melt (t(BAI), id.vars = id.vars,variable_name=Treename)




#####################################################
BAI10yr <- apply(X=BAI[,2:11],MARGIN=1,FUN=mean) 
BAI10yrtot <- apply (X=BAI[,2:11],MARGIN=1,FUN=sum)
BAI10yrmax <- apply (X=BAI[,2:11],MARGIN=1,FUN=max)
BAI10yrmin <- apply (X=BAI[,2:11],MARGIN=1,FUN=min)
tree_info$BAI10yr <-BAI10yr 
write.csv(tree_info, "PIPO-BAI10yr_v2.csv")#####################
######################################################

levels(tree_info$Band) <- list(L="L", M = "M", H = "H")
levels(tree_info$Comp) <- list(H = "H", L="L")


# --------- Preliminary BAI Boxplot --------------------
# getting Low, Med, High to show up right
#levels(tree_info$Band) <- list(L="L", M = "M", H = "H")
quartz(width=4, height=4)
par(mgp=c(2.4,1,0), mar=c(4,4,2,1), cex.lab=1.2)
boxplot(BAI10yrmin~Comp*Band, data=tree_info, col=c("grey", "white"),boxwex=.7, outwex=0
        , at = c(1,2,4,5,7,8)
        , ylab = expression("10yr Mean Basal Area Increment (m" * m^2 * ")")
        , xlab = "Elevation Band"
        , xaxt = "n"
        , main = "MIN"
        )
axis(1,at=c(1.5,4.5,7.5), labels=c("Low", "Mid", "High"))
legend("topleft",legend=c("Competitive", "Noncompetitive"),fill=c("grey", "white"), bty="n")

summary(aov(BAI10yr~Comp*Band, data=tree_info))
### Yay! the comp, elevation, and interaction effect is still significant with BAI


##########################################################
#                 more competition analysis     
#        analyzing differences and continuous metrics
##########################################################

#---calculating the difference in BAI between pairs
# first getting all of the pairs w/comp treatment
tree_info$Pair[which(tree_info$Pair==1)] <- "01" #that pesky 1 vs 01 coming up again
pairs <- unique(paste(tree_info$Species, tree_info$Band, tree_info$Pair, sep="-"))
BAI_dif <- c()
N_Cr_dif <- c()
BA_tot_dif <- c()
BA_totL_dif <- c()
BA_totLR_dif <- c()
BA_same_dif <- c()
BA_sameL_dif <- c()
BA_sameLR_dif <- c()
BA_other_dif <- c()
BA_otherL_dif <- c()
BA_otherLR_dif <- c()
DBH_dif <- c()
in5_tot_dif <- c()
in5_L_dif <- c()
in5_same_dif <- c()
in5_sameL_dif <- c()
in5_other_dif <- c()
in5_otherL_dif <- c()
for (i in 1:length(pairs)) {
  tmp <- tree_info[grep(pattern=pairs[i], x=tree_info$Tag),]
  #print(pairs[i]) # in for trouble shooting purposes
  #print(tmp$BAI10yr[tmp$Comp=="L"] - tmp$BAI10yr[tmp$Comp=="H"])
  BAI_dif[i] <- tmp$BAI10yr[tmp$Comp=="L"] - tmp$BAI10yr[tmp$Comp=="H"]
  N_Cr_dif[i] <- tmp$N_Cr[tmp$Comp=="L"] - tmp$N_Cr[tmp$Comp=="H"]
  BA_tot_dif[i] <- tmp$BA_tot[tmp$Comp=="L"] - tmp$BA_tot[tmp$Comp=="H"]
  BA_totL_dif[i] <- tmp$BA_totL[tmp$Comp=="L"] - tmp$BA_totL[tmp$Comp=="H"]
  BA_totLR_dif[i] <- (tmp$BA_totL[tmp$Comp=="L"] + tmp$BA_totR[tmp$Comp=="L"]) - (tmp$BA_totL[tmp$Comp=="H"] + tmp$BA_totR[tmp$Comp=="H"])
  BA_sameL_dif[i] <- tmp$BA_sameL[tmp$Comp=="L"] - tmp$BA_sameL[tmp$Comp=="H"]
  BA_sameLR_dif[i] <- tmp$BA_sameLR[tmp$Comp=="L"] - tmp$BA_sameLR[tmp$Comp=="H"]
  BA_other_dif[i] <- tmp$BA_other_all[tmp$Comp=="L"] - tmp$BA_other_all[tmp$Comp=="H"]
  BA_otherL_dif[i] <- tmp$BA_otherL[tmp$Comp=="L"] - tmp$BA_otherL[tmp$Comp=="H"]
  BA_otherLR_dif[i] <- tmp$BA_otherLR[tmp$Comp=="L"] - tmp$BA_otherLR[tmp$Comp=="H"]  
  DBH_dif[i] <- tmp$DBH[tmp$Comp=="L"] - tmp$DBH[tmp$Comp=="H"]
  in5_tot_dif[i] <- tmp$in5_tot[tmp$Comp=="L"] - tmp$in5_tot[tmp$Comp=="H"]
  in5_L_dif[i] <- tmp$in5_L[tmp$Comp=="L"] - tmp$in5_L[tmp$Comp=="H"]
  in5_same_dif[i] <- tmp$in5_same[tmp$Comp=="L"] - tmp$in5_same[tmp$Comp=="H"]
  in5_sameL_dif[i] <- tmp$in5_sameL[tmp$Comp=="L"] - tmp$in5_sameL[tmp$Comp=="H"]
  in5_other_dif[i] <- tmp$in5_other[tmp$Comp=="L"] - tmp$in5_other[tmp$Comp=="H"]
  in5_otherL_dif[i] <- tmp$in5_otherL[tmp$Comp=="L"] - tmp$in5_otherL[tmp$Comp=="H"]
  
}
# make elevation factor
elev <- c(rep("H", times=12), rep("L", times=10), rep("M", times=10))
difs <- data.frame(cbind(pairs, elev,BAI_dif,
              N_Cr_dif,
              BA_tot_dif,
              BA_totL_dif,
              BA_totLR_dif,
              BA_same_dif,
              BA_sameL_dif, BA_sameLR_dif,
              BA_other_dif, 
              BA_otherL_dif,
              BA_otherLR_dif,
              DBH_dif,
              in5_tot_dif,
              in5_L_dif,
              in5_same_dif,
              in5_sameL_dif,
              in5_other_dif,
              in5_otherL_dif))
ncol(difs)

# plotting all the variables at once
quartz(width=9, height=6.8)
par(mar=c(3,3,2, 1), mfrow=c(3,5))
cols <- brewer.pal(n=3, name="Set1")
for(i in 4:(ncol(difs)-1)){
  var <- colnames(difs)[i]
  plot(BAI_dif~as.numeric(as.character(difs[,i])), type="n", main=var)
  abline(h=0, lty=2)
  points(BAI_dif~as.numeric(as.character(difs[,i])), subset=which(difs$elev=="H"), pch=16, col=cols[1])
  points(BAI_dif~as.numeric(as.character(difs[,i])), subset=which(difs$elev=="M"), pch=16, col=cols[2])
  points(BAI_dif~as.numeric(as.character(difs[,i])), subset=which(difs$elev=="L"), pch=16, col=cols[3])
  H <- lm(BAI_dif~as.numeric(as.character(difs[,i])), subset=which(difs$elev=="H"))
  M <- lm(BAI_dif~as.numeric(as.character(difs[,i])), subset=which(difs$elev=="M"))
  L <- lm(BAI_dif~as.numeric(as.character(difs[,i])), subset=which(difs$elev=="L"))
  abline(H, lwd=2, col=cols[1])
  abline(M, lwd=2, col=cols[2])
  abline(L, lwd=2, col=cols[3])
  mtext(side=1, padj=-.35, cex=.6, text=paste("High p =", round(coef(summary(H))[2,4],digits=3),"\n", "Med p=", round(coef(summary(M))[2,4], digits=3), "\n", "Low p=", round(coef(summary(L))[2,4], digits=3)))

}


# plotting individual variables
quartz(width=6, height=4.5)
par(mfrow=c(1,3))
plot(BAI_dif~BA_totL_dif, type="n", main="total Live")
abline(h=0, lty=2, col="gray")
points(BAI_dif~BA_totL_dif, subset=which(elev=="H"), pch=16, col="blue")
points(BAI_dif~BA_totL_dif, subset=which(elev=="M"), pch=16, col="purple")
points(BAI_dif~BA_totL_dif, subset=which(elev=="L"), pch=16, col="red")
H <- lm(BAI_dif~BA_totL_dif, subset=which(elev=="H"))
M <- lm(BAI_dif~BA_totL_dif, subset=which(elev=="M"))
L <- lm(BAI_dif~BA_totL_dif, subset=which(elev=="L"))
abline(H, lwd=2, col="blue")
abline(M, lwd=2, col="purple")
abline(L, lwd=2, col="red")
cat("High p =", coef(summary(H))[2,4], "\n", "Med p=", coef(summary(M))[2,4], "\n", "Low p=", coef(summary(L))[2,4], "\n")








#===============================================================
############## POTR Analysis ###################################
#===============================================================


#++++++++++++++++++++++++++++++++++++++++++++++
# NOTE: code between the ++++ markes resuses
# nomenclature for different species. so must be rerun for each species
## low elevation
csv1 <- "CO-POTR-L_RWdata01_29_14.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
treedata1 <- core_average(data=data1)
## mid elevation
csv2 <- "CO-POTR-M_RWdata01_31_14.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
treedata2 <- core_average(data=data2)
## Somehow some of my cores are labeled "PIPO" in this spreadsheet.
treedata2$Species <- rep("POTR", times=nrow(treedata2))

## High elevation
csv3 <- "CO-POTR-H_RWdata04_14_14v2.csv" # 6/22/15 had to get rid of a label H-02L
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3) 

#------------- Getting 10yr Average Growth ---------------
# tree10yr <- rbind(treedata1[,1:15], treedata2[1:15], treedata3[,1:15])
# tenyrmean <- apply (X=tree10yr[,5:15],MARGIN=1, FUN=mean) # getting 10 yr average growth
# Treat <- rep(c("C", "N"),times= length(tenyrmean)/2) 
# #can't think of a better way to get competitive treatment right now
# # but seems to work because H is always first. switching to C and N to avoid confusion
# 
# comp_POTR <- cbind(tree10yr[1:4], Treat, tenyrmean)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# test <- aov (tenyrmean~Treat*Elev, data = comp_POTR)
# TukeyHSD(test)
# 
# # --------- Preliminary Boxplot --------------------
# quartz(width=4, height=4.5)
# boxplot(tenyrmean~Treat*Elev, data=comp_POTR, col=c("grey", "white")
#         , ylab = "Average Ring Width (mm)"
#         , xlab = "Elevation"
#         , xaxt = "n"
#         , main = "POTR 10yr average growth"
#         , at=c(1,2,4,5,7,8))
# axis(1,at=c(1.5,4.5,7.5), labels=c("Low", "Mid", "High"))
# legend("topleft",legend=c("Competitive", "Noncompetitive"),fill=c("grey", "white"), bty="n")
# 
# # interaction plot
# interaction.plot(x.factor=comp_PIPO$Elev, trace.factor=comp_PIPO$Treat, response=comp_PIPO$tenyrmean)
# 


# make tags for common reference
ncols <- min(c(ncol(treedata1), ncol(treedata2), ncol(treedata3)))
rwdata_all <- rbind(treedata1[,1:ncols], treedata2[1:ncols], treedata3[1:ncols])
rwdata_all$Tree <- as.character(rwdata_all$Tree)
####### Need to change when switching nomenclature
# rwdata_all$Tree[which(rwdata_all$Tree=="1L")] <- "01L" #having problems with 1 vs 10,11,12 later in script, so changing to 01
# rwdata_all$Tree[which(rwdata_all$Tree=="1H")] <- "01H"
rwdata_all$Tag <- paste(rwdata_all$Species, rwdata_all$Elev, rwdata_all$Tree, sep = "-")
#made unique tree tag
rwdata_alltag <- rwdata_all[,c(1:4,ncol(rwdata_all),5:(ncol(rwdata_all)-1))]
#reordered columns to get tag in right place
rwdata_alltag <- rwdata_alltag[order(rwdata_alltag$Tag),]

# make a treeDBH data frame with just the tag, DBH and bark
tree_info <- COtreesnew[which(COtreesnew$Tag %in% rwdata_alltag$Tag),]
#tree_info <- COtreesnew[match(rwdata)]
#MtreeDBH <- COtreesnew[which(COtreesnew$Tag %in% treedata2$Tag), c(6,12,13)]
#LtreeDBH <- COtreesnew[which(COtreesnew$Tag %in% treedata1$Tag), c(6,12,13)]
tree_info <- tree_info[order(tree_info$Tag),] # get in in the same order as treedata

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


tree_infoPOTR <- tree_info
rwdata_alltagPOTR <- rwdata_alltag

rwdata_alltagPOTR$Tag==tree_infoPOTR$Tag


#+++++++++++++++++++++++++++++++++++++++
BAI <- basal_area_inc(rwdata_alltagPOTR, tree_infoPOTR)
min(BAI, na.rm=T) # Yeesh! I've got some trees that are off by quite a ways. worst is off by 1.5cm



#####################################################
      # averaging 2012-2003
BAI10yr <- apply(X=BAI[,2:11],MARGIN=1,FUN=mean) 
BAI10yrtot <- apply (X=BAI[,2:11],MARGIN=1,FUN=sum)
BAI10yrmax <- apply (X=BAI[,2:11],MARGIN=1,FUN=max)
BAI10yrmin <- apply (X=BAI[,2:11],MARGIN=1,FUN=min)
tree_infoPOTR$BAI10yr <-BAI10yr 
write.csv(tree_infoPOTR, "POTR-BAI10yr_v1.csv")#####################
######################################################








#===============================================================
############## ABLA Analysis ###################################
#===============================================================
csv1 <- "CO-ABLA-L-RWdata05_02_14_v2.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
colnames(data1)[1] <- "Transect" 
#coredata1 <- core_clean(data=data1) #unaveraged cores (i.e. 2 cores per tree still)
treedata1 <- core_average(data=data1) # cores averaged per tree
## mid elevation
csv2 <- "CO-ABLA-M-RWdata05_05_14.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
colnames(data2)[1] <- "Transect" 
#coredata2 <- core_clean(data=data2)
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "CO-ABLA-H-RWdata04_30_14.csv"
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
#coredata3 <- core_clean(data=data3)
treedata3 <- core_average(data=data3)

# make tags for common reference
ncols <- min(c(ncol(treedata1), ncol(treedata2), ncol(treedata3)))
rwdata_all <- rbind(treedata1[,1:ncols], treedata2[1:ncols], treedata3[1:ncols])
rwdata_all$Tree <- as.character(rwdata_all$Tree)
####### Need to change when switching nomenclature
# rwdata_all$Tree[which(rwdata_all$Tree=="1L")] <- "01L" #having problems with 1 vs 10,11,12 later in script, so changing to 01
# rwdata_all$Tree[which(rwdata_all$Tree=="1H")] <- "01H"
rwdata_all$Tag <- paste(rwdata_all$Species, rwdata_all$Elev, rwdata_all$Tree, sep = "-")
#made unique tree tag
rwdata_alltag <- rwdata_all[,c(1:4,ncol(rwdata_all),5:(ncol(rwdata_all)-1))]
#reordered columns to get tag in right place
rwdata_alltag <- rwdata_alltag[order(rwdata_alltag$Tag),]

# make a treeDBH data frame with just the tag, DBH and bark
tree_info <- COtreesnew[which(COtreesnew$Tag %in% rwdata_alltag$Tag),]
#tree_info <- COtreesnew[match(rwdata)]
#MtreeDBH <- COtreesnew[which(COtreesnew$Tag %in% treedata2$Tag), c(6,12,13)]
#LtreeDBH <- COtreesnew[which(COtreesnew$Tag %in% treedata1$Tag), c(6,12,13)]
tree_info <- tree_info[order(tree_info$Tag),] # get in in the same order as treedata

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


tree_infoABLA <- tree_info
rwdata_alltagABLA <- rwdata_alltag

rwdata_alltagABLA$Tag==tree_infoABLA$Tag



############################# Save my PIPO and my POTR rws and tree infos ####################
# write.csv(rwdata_alltagPIPO, "PIPO_avgRingWidth_2-19-15.csv")
# write.csv(tree_infoPIPO, "PIPO_TreeInfo.csv")
# write.csv(rwdata_alltagPOTR, "POTR_avgRingWidth_2-19-15.csv")
# write.csv(tree_infoPOTR, "POTR_TreeInfo.csv")
# write.csv(rwdata_alltagABLA, "ABLA_avgRingWidth_2-19-15.csv")
# write.csv(tree_infoABLA, "ABLA_TreeInfo.csv")
###############################################################################################




