######################################################
#       Script for producing density plots 
#       of pairwise correlation coefficients
#            for all three species in CO
#######################################################


# because I can't think of any easy way to get all the 
# data in without overwriting things i've stored in objects of the same name
# I've just copied all the requisite code for each species. It's ugly, but oh well


####### needed functions
source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/fn_core_trunc.R")
source("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/CO_RingWidths/CO_RW_csvs_11-13/fn_core_average.R")


####### making plot
xlim <- c(-0.5,1)
ylim <- c(0, 4.7)

quartz(width=5, heigh=6.5)
par(mfrow=c(3,1), mar=c(3,5,0,0), oma=c(2,0,1,1))

### Starting with ABLA on top
#---Getting DATA------------------------
## low elevation
csv1 <- "CO-ABLA-L-RWdata05_02_14_v2.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
colnames(data1)[1] <- "Transect" 
treedata1 <- core_average(data=data1) # cores averaged per tree
## mid elevation
csv2 <- "CO-ABLA-M-RWdata05_05_14.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
colnames(data2)[1] <- "Transect" 
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "CO-ABLA-H-RWdata04_30_14.csv"
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3)

### customize these to get dates
# for PRISM data from climate WNA, first full water year = 1902, last year =2011
#workingdata <-treedata1
start <- "1902"   # first date for time series
end <- "2011"     # last date of time series
method <- "Spline" # detrending method desired. c("Spline", "ModNegExp", "Mean")
dropyrs <- 3
####

# High Elevation
rw.detH <- core_trunc(workingdata=treedata3, detmethod=method, detrend=TRUE, drop=dropyrs)
# Mid Elevation
rw.detM <- core_trunc(workingdata=treedata2, detmethod=method, detrend=TRUE, drop=dropyrs)
# Low Elevation
rw.detL <- core_trunc(workingdata=treedata1, detmethod=method, detrend=TRUE, drop=dropyrs)

tmpdata <- rw.detH # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.H <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.H<-c(mean(z.H[,"cor"], na.rm=TRUE),sd(z.H[,"cor"], na.rm=TRUE), ncol(tmpdata))

###---- Mid elevation
tmpdata <- rw.detM # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.M <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.M<-c(mean(z.M[,"cor"], na.rm=TRUE),sd(z.M[,"cor"], na.rm=TRUE), ncol(tmpdata))

###--- Low elevation
tmpdata <- rw.detL # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.L <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.L<-c(mean(z.L[,"cor"], na.rm=TRUE),sd(z.L[,"cor"], na.rm=TRUE), ncol(tmpdata))

corr.all <- rbind(z.L,z.M,z.H)
elevs <- c(rep("L", times=nrow(z.L)), rep("M", times=nrow(z.M)), rep("H", times=nrow(z.H)))
corr.all <- cbind(corr.all, elevs)
#-------------------Plotting ABLA-----------------

plot(density(corr.all[,"cor"]), type="n",xlim=xlim, ylim =ylim, main = "", xlab="", yaxs = "i", ylab="")
abline(v=tapply(X=corr.all$cor,INDEX=elevs,FUN=mean), lty=2, lwd=3, col=c('blue3',"red3", "purple3"))

polygon(x=c(density(z.H$cor)$x,min(density(z.H$cor)$x)), y=c(density(z.H$cor)$y,0), border="black", col =rgb(0,0,205,150,max=255))
polygon(x=c(density(z.M$cor)$x,min(density(z.M$cor)$x)), y=c(density(z.M$cor)$y,0), border="black", col= rgb(125,38,205,180, max=255))
polygon(x=c(density(z.L$cor)$x,min(density(z.L$cor)$x)), y=c(density(z.L$cor)$y,0), border="black", col=rgb(205,0,0,150, max=255) )
legend("left",legend=c("High", "Mid", "Low"),fill=c(rgb(0,0,205,150,max=255), rgb(125,38,205,180, max=255),rgb(205,0,0,150, max=255)), border="black",title="Elevation", bty="n")
text(x=-0.5, y=ylim[2]-.5,labels="Abies lasiocarpa",font=4, cex=1.5, pos=4)
#abline(v=tapply(X=corr.all$cor,INDEX=elevs,FUN=mean), lty=2, lwd=3, col=c('blue3',"red3", "purple3"))

### Heading down hill to POTR
#---------Getting DATA-------
csv1 <- "CO-POTR-L_RWdata01_29_14.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
treedata1 <- core_average(data=data1) # cores averaged per tree
## mid elevation
csv2 <- "CO-POTR-M_RWdata01_31_14.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "CO-POTR-H_RWdata04_14_14v2.csv"
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3)

# High Elevation
rw.detH <- core_trunc(workingdata=treedata3, detmethod=method, detrend=TRUE, drop=dropyrs)
# Mid Elevation
rw.detM <- core_trunc(workingdata=treedata2, detmethod=method, detrend=TRUE, drop=dropyrs)
# Low Elevation
rw.detL <- core_trunc(workingdata=treedata1, detmethod=method, detrend=TRUE, drop=dropyrs)

tmpdata <- rw.detH # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.H <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.H<-c(mean(z.H[,"cor"], na.rm=TRUE),sd(z.H[,"cor"], na.rm=TRUE), ncol(tmpdata))

###---- Mid elevation
tmpdata <- rw.detM # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.M <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.M<-c(mean(z.M[,"cor"], na.rm=TRUE),sd(z.M[,"cor"], na.rm=TRUE), ncol(tmpdata))

###--- Low elevation
tmpdata <- rw.detL # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.L <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.L<-c(mean(z.L[,"cor"], na.rm=TRUE),sd(z.L[,"cor"], na.rm=TRUE), ncol(tmpdata))

# one big df for all the pairwise stuff
corr.all <- rbind(z.L,z.M,z.H)
elevs <- c(rep("L", times=nrow(z.L)), rep("M", times=nrow(z.M)), rep("H", times=nrow(z.H)))
corr.all <- cbind(corr.all, elevs)


#--------- Plotting POTR ----------------
plot(density(corr.all[,"cor"]), type="n",xlim=xlim, ylim =ylim, main = "", xlab="", yaxs = "i")
abline(v=tapply(X=corr.all$cor,INDEX=elevs,FUN=mean), lty=2, lwd=3, col=c('blue3',"red3", "purple3"))
polygon(x=c(density(z.H$cor)$x,min(density(z.H$cor)$x)), y=c(density(z.H$cor)$y,0), border="black", col =rgb(0,0,205,150,max=255))
polygon(x=c(density(z.M$cor)$x,min(density(z.M$cor)$x)), y=c(density(z.M$cor)$y,0), border="black", col= rgb(125,38,205,180, max=255))
polygon(x=c(density(z.L$cor)$x,min(density(z.L$cor)$x)), y=c(density(z.L$cor)$y,0), border="black", col=rgb(205,0,0,150, max=255) )
text(x=-0.5, y=ylim[2]-.5,labels="Populus tremuloides",font=4, cex=1.5, pos=4)
#abline(v=tapply(X=corr.all$cor,INDEX=elevs,FUN=mean), lty=2, lwd=3, col=c('blue3',"red3", "purple3"))


### And then donw all the way to PIPO
#----------Getting DATA ---------------
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

# High Elevation
rw.detH <- core_trunc(workingdata=treedata3, detmethod=method, detrend=TRUE, drop=dropyrs)
# Mid Elevation
rw.detM <- core_trunc(workingdata=treedata2, detmethod=method, detrend=TRUE, drop=dropyrs)
# Low Elevation
rw.detL <- core_trunc(workingdata=treedata1, detmethod=method, detrend=TRUE, drop=dropyrs)

tmpdata <- rw.detH # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.H <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.H<-c(mean(z.H[,"cor"], na.rm=TRUE),sd(z.H[,"cor"], na.rm=TRUE), ncol(tmpdata))

###---- Mid elevation
tmpdata <- rw.detM # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.M <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.M<-c(mean(z.M[,"cor"], na.rm=TRUE),sd(z.M[,"cor"], na.rm=TRUE), ncol(tmpdata))

###--- Low elevation
tmpdata <- rw.detL # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.L <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.L<-c(mean(z.L[,"cor"], na.rm=TRUE),sd(z.L[,"cor"], na.rm=TRUE), ncol(tmpdata))

corr.all <- rbind(z.L,z.M,z.H)
elevs <- c(rep("L", times=nrow(z.L)), rep("M", times=nrow(z.M)), rep("H", times=nrow(z.H)))
corr.all <- cbind(corr.all, elevs)

#---------- Plotting PIPO -----------------

plot(density(corr.all[,"cor"]), type="n",xlim=xlim, ylim =ylim, main = "", xlab="", yaxs = "i", ylab="")
abline(v=tapply(X=corr.all$cor,INDEX=elevs,FUN=mean), lty=2, lwd=3, col=c('blue3',"red3", "purple3"))
polygon(x=c(density(z.H$cor)$x,min(density(z.H$cor)$x)), y=c(density(z.H$cor)$y,0), border="black", col =rgb(0,0,205,150,max=255))
polygon(x=c(density(z.M$cor)$x,min(density(z.M$cor)$x)), y=c(density(z.M$cor)$y,0), border="black", col= rgb(125,38,205,180, max=255))
polygon(x=c(density(z.L$cor)$x,min(density(z.L$cor)$x)), y=c(density(z.L$cor)$y,0), border="black", col=rgb(205,0,0,150, max=255) )
text(x=-0.5, y=ylim[2]-.5,labels="Pinus ponderosa",font=4, cex=1.5, pos=4)
#abline(v=tapply(X=corr.all$cor,INDEX=elevs,FUN=mean), lty=2, lwd=3, col=c('blue3',"red3", "purple3"))

mtext(text="Tree Pairwise Pearson Correlation Coefficients", side=1,outer=T, font=1)






########################################################
########   Trying boxplot option #######################
########################################################


#-------------- Getting ABLA data -------------
## low elevation
csv1 <- "CO-ABLA-L-RWdata05_02_14_v2.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
colnames(data1)[1] <- "Transect" 
treedata1 <- core_average(data=data1) # cores averaged per tree
## mid elevation
csv2 <- "CO-ABLA-M-RWdata05_05_14.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
colnames(data2)[1] <- "Transect" 
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "CO-ABLA-H-RWdata04_30_14.csv"
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3)

### customize these to get dates
# for PRISM data from climate WNA, first full water year = 1902, last year =2011
#workingdata <-treedata1
start <- "1902"   # first date for time series
end <- "2011"     # last date of time series
method <- "Spline" # detrending method desired. c("Spline", "ModNegExp", "Mean")
dropyrs <- 3
####

# High Elevation
rw.detH <- core_trunc(workingdata=treedata3, detmethod=method, detrend=TRUE, drop=dropyrs)
# Mid Elevation
rw.detM <- core_trunc(workingdata=treedata2, detmethod=method, detrend=TRUE, drop=dropyrs)
# Low Elevation
rw.detL <- core_trunc(workingdata=treedata1, detmethod=method, detrend=TRUE, drop=dropyrs)

tmpdata <- rw.detH # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.H <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.H<-c(mean(z.H[,"cor"], na.rm=TRUE),sd(z.H[,"cor"], na.rm=TRUE), ncol(tmpdata))

###---- Mid elevation
tmpdata <- rw.detM # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.M <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.M<-c(mean(z.M[,"cor"], na.rm=TRUE),sd(z.M[,"cor"], na.rm=TRUE), ncol(tmpdata))

###--- Low elevation
tmpdata <- rw.detL # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.L <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.L<-c(mean(z.L[,"cor"], na.rm=TRUE),sd(z.L[,"cor"], na.rm=TRUE), ncol(tmpdata))

corr.all <- rbind(z.L,z.M,z.H)
elevs <- c(rep("L", times=nrow(z.L)), rep("M", times=nrow(z.M)), rep("H", times=nrow(z.H)))
species <- rep("ABLA", times=length(elevs))
corr.all.ABLA <- cbind(corr.all, elevs, species)


#---------Getting POTR DATA-------
csv1 <- "CO-POTR-L_RWdata01_29_14.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
treedata1 <- core_average(data=data1) # cores averaged per tree
## mid elevation
csv2 <- "CO-POTR-M_RWdata01_31_14.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "CO-POTR-H_RWdata04_14_14v2.csv"
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3)

# High Elevation
rw.detH <- core_trunc(workingdata=treedata3, detmethod=method, detrend=TRUE, drop=dropyrs)
# Mid Elevation
rw.detM <- core_trunc(workingdata=treedata2, detmethod=method, detrend=TRUE, drop=dropyrs)
# Low Elevation
rw.detL <- core_trunc(workingdata=treedata1, detmethod=method, detrend=TRUE, drop=dropyrs)

tmpdata <- rw.detH # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.H <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.H<-c(mean(z.H[,"cor"], na.rm=TRUE),sd(z.H[,"cor"], na.rm=TRUE), ncol(tmpdata))

###---- Mid elevation
tmpdata <- rw.detM # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.M <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.M<-c(mean(z.M[,"cor"], na.rm=TRUE),sd(z.M[,"cor"], na.rm=TRUE), ncol(tmpdata))

###--- Low elevation
tmpdata <- rw.detL # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.L <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.L<-c(mean(z.L[,"cor"], na.rm=TRUE),sd(z.L[,"cor"], na.rm=TRUE), ncol(tmpdata))

# one big df for all the pairwise stuff
corr.all <- rbind(z.L,z.M,z.H)
elevs <- c(rep("L", times=nrow(z.L)), rep("M", times=nrow(z.M)), rep("H", times=nrow(z.H)))
species <- rep("POTR", times=length(elevs))
corr.all.POTR <- cbind(corr.all, elevs, species)

#----------Getting PIPO DATA ---------------
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

# High Elevation
rw.detH <- core_trunc(workingdata=treedata3, detmethod=method, detrend=TRUE, drop=dropyrs)
# Mid Elevation
rw.detM <- core_trunc(workingdata=treedata2, detmethod=method, detrend=TRUE, drop=dropyrs)
# Low Elevation
rw.detL <- core_trunc(workingdata=treedata1, detmethod=method, detrend=TRUE, drop=dropyrs)

tmpdata <- rw.detH # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.H <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.H<-c(mean(z.H[,"cor"], na.rm=TRUE),sd(z.H[,"cor"], na.rm=TRUE), ncol(tmpdata))

###---- Mid elevation
tmpdata <- rw.detM # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.M <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.M<-c(mean(z.M[,"cor"], na.rm=TRUE),sd(z.M[,"cor"], na.rm=TRUE), ncol(tmpdata))

###--- Low elevation
tmpdata <- rw.detL # designate which data to use
# get pairwise correlation matrix
r<- cor(tmpdata, use = "na.or.complete") # get corr matrix
y<- which(lower.tri(r), TRUE) # get values to subset to lower triangular
z.L <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y]) # make a dataframe of each pairwise correlation
# make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
corr.L<-c(mean(z.L[,"cor"], na.rm=TRUE),sd(z.L[,"cor"], na.rm=TRUE), ncol(tmpdata))

corr.all <- rbind(z.L,z.M,z.H)
elevs <- c(rep("L", times=nrow(z.L)), rep("M", times=nrow(z.M)), rep("H", times=nrow(z.H)))
species <- rep("PIPO", times=length(elevs))
corr.all.PIPO <- cbind(corr.all, elevs,species)

#----------- plotting boxplot ----------------------
correlations <- rbind(corr.all.ABLA,corr.all.POTR,corr.all.PIPO)
levels(correlations$elevs) <- list(L = "L", M = "M", H="H")
quartz(width=4, height=6)
par(mar=c(4,2,2,1))
boxplot(cor~elevs,data=correlations,subset=species=="PIPO",horizontal=T, xlim=c(0,18), ylim=c(-0.5, 1), col=c("red3", "purple3","blue3"), yaxt="n")
boxplot(cor~elevs,data=correlations,subset=species=="POTR",horizontal=T, add=T, at=c(7,8,9), col=c("red3", "purple3","blue3"), yaxt="n")
boxplot(cor~elevs,data=correlations,subset=species=="ABLA",horizontal=T, add=T, at=c(13,14,15), col=c("red3", "purple3","blue3"), yaxt="n")
text(x=c(-.5,-.5,-.5), y=c(5,11,17),labels=c("Pinus ponderosa", "Populus tremuloides", "Abies lasiocarpa"), font=4, pos=4)
mtext(text="Tree Pairwise Pearson Correlation Coefficients", side=1,line=2.5, font=1)
legend("bottomleft",legend=c("High", "Mid", "Low"),fill=c("blue3","purple3","red3"), border="black",title="Elevation", bty="n", cex=.9)#c(rgb(0,0,205,150,max=255), rgb(125,38,205,180, max=255),rgb(205,0,0,150, max=255)), border="black",title="Elevation", bty="n")



summary(aov(cor~elevs,data=correlations,subset=species=="ABLA")) # p=0.616
summary(aov(cor~elevs,data=correlations,subset=species=="POTR")) # p= <2e-16
summary(aov(cor~elevs,data=correlations,subset=species=="PIPO")) # p <2e-16

TukeyHSD(aov(cor~elevs,data=correlations,subset=species=="POTR")) # so they're all different
