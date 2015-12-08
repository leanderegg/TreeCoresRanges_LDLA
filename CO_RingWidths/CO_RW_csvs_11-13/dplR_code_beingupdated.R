####################################################
########### Crossdating with dplR ##################
###################################################

# set wd to source location
   #NOTE: this will need to be redirected to /MT_RingWidths... when we move to the next state

 setwd("C:/Users/Western Hemlock/Dropbox/Postdocs&GradProjects/TreeCoresRanges(Leander)/MT_RingWidths/MT_RW_csvs_05-14")
 setwd("C:/Users/Western Hemlock/Dropbox/Postdocs&GradProjects/TreeCoresRanges(Leander)/WA_RingWidths/WA_RW_csvs")
 setwd("C:/Users/Western Hemlock/Dropbox/Postdocs&GradProjects/TreeCoresRanges(Leander)/CO_RingWidths/CO_RW_csvs_11-13")


##### RUN this Code at beginning ######################################
# loads the R package we need
library(dplR)




#----------------------------------------------------------------------------
#--------------------- Crossdate function to create plot --------------
#--------------------------------------------------------------
#----------------------------------------------------
crossdate <- function(data)
{
data$Analysis.Date.Time <- as.POSIXct(strptime(data$Analysis.Date.Time, format = "%m/%d/%Y %H:%M"))
un <- unique(data$TreeID)
dataclean <- data[1:length(un),]
for(i in 1:length(un)){
data_subset <- data[which(data$TreeID==un[i]),]
print(paste(un[i], nrow(data_subset[which(data_subset$Analysis.Date.Time  == max(data_subset$Analysis.Date.Time)),])))
dataclean[i,] <- data_subset[which(data_subset$Analysis.Date.Time  == max(data_subset$Analysis.Date.Time)),]
  }
#dataclean$Tag

coredataraw <- dataclean




badcols <- c( "Path.identification","Site.identification" ,"YearLastRing", "Sapwood","Tree.height","Tree.age"
	,"SectionHeight" , "User.variable" ,"DataType","OffsetToNext","ImageName"
	, "Acquisition.Date.Time","Modified.Date.Time", "ImageSizeH.V.NBits.Channel"
	,"CalibMethod.XCal.YCal.EditedDendro","ScannerCamera.Make.Model.Software","LensFocLength..35mm."
	,"PathBegX.BegY.EndX.EndY.Width","RingBoundary.AutoMan.Meth.Precise","EarlywoodDef"
	,"DensActive.Media.Calib","DensNSteps.MatDens.Interpol.SampleThick","DensStepsThick"
	,"DensStepsLightInt","DensStepsWoodDens","DiskArea","DiskPerim","DiskAvgDiam"
	,"DiskFormCoef","CompWoodArea","VoidArea","PathLength")
 # columns to kill from Windendro auto output
coredata <- coredataraw[,!(colnames(coredataraw) %in% badcols)]
coredata[,"X"] <- NULL # removing 2013 (year core was collected)


	spstdat<-coredata
     		spstdat.rwl<-t(spstdat[,10:ncol(spstdat)]) ### replace with second year of measurements (deleting 2013 and starting with 2012)
     		rownames(spstdat.rwl)<-c(2012:(2012-nrow(spstdat.rwl)+1))
     		colnames(spstdat.rwl)<-spstdat$TreeID
     		spstdat.rwl.df<-as.data.frame(spstdat.rwl)
     	spst.stats<-corr.rwl.seg(spstdat.rwl.df,seg.length=20, bin.floor=round(min(as.numeric(rownames(spstdat.rwl))),digits=-1)+5, main = data$Elev[1])
#seg.length can be changed to analyze cores on a finer scale (right now set at 20 years)		
return(spstdat.rwl.df)
}

#-----------------------------
#-------------------------------------------
#------------------------------------------------------------------
#################### everything above here, just run once upon startup








#==============================================================
#---------------------- Analysis ------------------------------
#==============================================================



#=============== Cleanding .csv up ============
# Note: before bringing into R, need to add: Transect, Species
		#, Elev, Tree, Core, Tag and TreeID columns in csv

#### cleaning the csv from WinDendro ####
# insert requisite columns and delete first row of csv
#csv1 <- "CO-PIPO-L_RWdata11_7_13test2.csv"
#csv2 <- "CO-PIPO-M_RWdata11_7_13.csv"
#csv3 <- "CO-PIPO-H_RWdata10_31_13.csv"
#csv <- "CO-PIPO-RWcomb.csv"

csv1 <- "CO-POTR-L_RWdata01_29_14.csv"
csv2 <- "CO-POTR-M_RWdata01_31_14.csv" 
csv3 <- "CO-POTR-H_RWdata04_14_14v2.csv"      

#csv1 <- "CO-ABLA-L-RWdata05_02_14_v2.csv"
#csv2 <- "CO-ABLA-M-RWdata05_05_14.csv"
#csv3 <- "CO-ABLA-H-RWdata04_30_14.csv"

#csv1 <- "MT-TSHE-L-RWdata05_07_14.csv"
#csv2 <- "MT-TSHE-M-RWdata05_12_14.csv"
#csv3 <- "MT-TSHE-H-RWdata05_09_14.csv"
 
#csv1 <- "MT-PSME-L-RWdata05_19_14.csv"
#csv2 <- "MT-PSME-M-RWdata06_04_14.csv"
#csv3 <- "MT-PSME-H-RWdata05_16_14.csv"

#csv1 <- "MT-ABLA-L-RWdata06_05_14.csv"
#csv2 <- "MT-ABLA-M-RWdata06_05_14.csv"
#csv3 <- "MT-ABLA-H-RWdata06_05_14.csv"
csv <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths//MT_RW_csvs_05-14/MT-ABLA-H-RWdata06_05_14.csv"

#csv2 <- "WA-TSHE-M-RWdata06_09_14.csv"
#csv3 <- "WA-TSHE-H-RWdata06_09_14.csv"

#csv1 <- "WA-PSME-L-RWdata06_09_14.csv"
#csv2 <- "WA-PSME-M-RWdata06_09_14.csv"
#csv3 <- "WA-PSME-H-RWdata06_10_14.csv"

#csv1 <- "WA-ABLA-L-RWdata06_10_14.csv"
#csv2 <- "WA-ABLA-M-RWdata06_10_14.csv"
#csv3 <- "WA-ABLA-H-RWdata06_13_14.csv"

#csv4 <- "CO-PIPO-Traits-ALL04_15_15.csv"
#csv5 <- "CO-PIPO-Traits-M-04_15_15.CSV"
#csv6 <- "CO-PIPO-Traits-H-04_15_15.CSV"

csv <- "CO-POTR-Traits-L-05_12_15.csv"
#csv <- "CO-POTR-Traits-M-05_08_15.csv"
#######################################
## Input file name
csv <- "CO-POTR-Lall.csv"    ### Change when reading in a new CSV
#########################################
#dataWP <- data1
#dataWA <- data3
 data1 <- read.csv(csv1, header=TRUE)
 data2 <- read.csv(csv2, header=TRUE)
 data3 <- read.csv(csv3, header=TRUE)
#names(data3)[1]<- "Transect"
#data4 <- read.csv(csv4, header=TRUE)
#data5 <- read.csv(csv5, header=TRUE)
#data6 <- read.csv(csv6, header=TRUE)

#head(data)
mincol <- 200


# data <- rbind(data1[1:mincol],data2[,1:mincol], data3[,1:mincol])#, data3[,1:mincol])
 data <- rbind(data3[1:mincol],data2[,1:mincol])#, data3[,1:mincol])
# data<- rbind(data4[,1:mincol], data5[,1:mincol],data6[1:mincol])
 data$TreeID <- with(data,paste(Elev,Tree,Core,sep="-"))

# reads csv into object "data"
# rerun this if you've just edited the csv file with new crossdating data from WinDendro
data <- read.csv(csv, header=TRUE)
# data$X <- NULL
# colnames(data)[1:mincol] <- colnames(data1)[1:mincol]
# data <- rbind(data[1:mincol],data1[,1:mincol])
#============= code for creating master chronology================
    # this line creates the plot with the red and blue horizontal lines showing all the cores
    # if you need to rerun this bit once you
spstdat.rwl.df<- crossdate(data=data)

#------------------------------------------------------
##===== pulling out individual cores to work with
#------------------------------------------------------

# isolate core of interest
trouble <- "L-E-T3-a"

# code that pulls out the individual core to analyze against the rest of the data set
k <- which(colnames(spstdat.rwl.df) == trouble) #input core of interest
flagtree<-colnames(spstdat.rwl.df)[k]
flagged=spstdat.rwl.df[,k]
names(flagged)=rownames(spstdat.rwl.df)
spstdat.rwl.df[,k]=NULL


# shows line graph of where correlation breaks down
seg.20<-corr.series.seg(rwl=spstdat.rwl.df,series=flagged,seg.length=20,bin.floor=05, main=paste(flagtree))#you may want to look at these line plots, but they use up too many windows if you are doing it for a whole stand
# shows ball and stick diagram that highlights additional or missing rings
ccf20=ccf.series.rwl(rwl=spstdat.rwl.df,series=flagged,seg.length=20,bin.floor=10,main=paste(trouble))
  

# SHIFT TO THE LEFT INDICATES EXTRA RING, SHIFT TO THE RIGHT INDICATES MISSING RING


data1 <- read.csv(csv1, header=TRUE)
data2 <- read.csv(csv2, header=TRUE)
data3 <- read.csv(csv3, header=TRUE)

Low<- crossdate(data=data1)
Mid <- crossdate(data=data2)
High <- crossdate(data=data3)

Ldet <- detrend(Low, method="Spline")
Lmaster <- chron(Ldet)

Mdet <- detrend(Mid, method="Spline")
Mmaster <- chron(Mdet)

Hdet <- detrend(High, method="Spline")
Hmaster <- chron(Hdet)

crn.plot(Lmaster)
crn.plot(Mmaster)
crn.plot(Hmaster)
lines(Mmaster[,1]~rownames(Mmaster), col="blue")
years <- as.numeric(rownames(Lmaster))
title <- "WA-PSME"
plot(Lmaster$xxxstd~years, col="darkred", type="l", xlim=c(1890,2010), lwd=2, main=title, ylim=c(0,1.5))
lines(Mmaster$xxxstd~rownames(Mmaster), col="darkgreen", lwd=2)
lines(Hmaster$xxxstd~rownames(Hmaster), col="blue", lwd=2)
grid()

newyears <- years-1
grid()
abline(v=2002)

trouble <- Low[,2]
series.rwl.plot(rwl=Low, series=trouble, series.yrs = as.numeric(row.names(Low)), seg.length=20, bin.floor=10)
write.crn(Lmaster,fname="TSME-L-Mastertest.crn")
write.crn(Mmaster, fname="TSME-M-Mastertest.crn")
write.rwl(rwl.df = Lmaster, fname="TSME-L-Mastertest.txt", format="tucson")
write.csv(Lmaster, "TSME-L-Master.csv")

### filling internal NAs for detrending
Midtest <- fill.internal.NA(Mid, fill="Spline")

## figuring out which cores have internal NAs
for (i in 1:ncol(Mid)){
  print(colnames(Mid)[i])
  tmp <-  detrend(data.frame(Mid[,i]), method="Spline")
}



######### ignore this code ##############

#=======================================================
## 			detrending 
#=======================================================


# looking at the trends for each core:
for(i in 1:ncol(spstdat.rwl.df)){windows()
detrend(rwl=spstdat.rwl.df[,c(i,i+1)], method="Spline", make.plot=T)}


## Creating a detrended rw data frame
rw.det <- detrend(rwl=spstdat.rwl.df, method="Spline")



### trying to just read in text file from WINDENDRO
#test <- read.table(file.choose(), header = TRUE, skip = 1, sep = "  "





# failing to get skeleton plot to work
skel.plot(spstdat.rwl.df[,2:3],yr.vec = as.numeric(rownames(spstdat.rwl.df)))
dim(spstdat.rwl.df)
head(spstdat.rwl.df)


test <- read.csv(file.choose(), header=FALSE, col.names=1)
head(test)




ccf20=ccf.series.rwl(rwl=spstdat.rwl.df,series=flagged,seg.length=20,bin.floor=0,main=paste(colnames(spstdat.rwl.df)[k]))





### Playing around with dplR functions

## the red noise power spectrum
library(graphics)
library(stats)
set.seed(123)
nyrs <- 500
yrs <- 1:nyrs

# Here is an ar1 time series with a mean of 2mm,
# an ar1 of phi, and sd of sigma
phi <- 0.7
sigma <- 0.3
sigma0 <- sqrt((1 - phi^2) * sigma^2)
x <- arima.sim(list(ar = phi), n = nyrs, sd = sigma0) + 2

# Here is a sine wave at f=0.1 to add in with an amplitude
# equal to half the sd of the red noise background
per <- 17 # change this to play with the period of the added sine wave
amp <- sigma0 / 2
wav <- amp * sin(2 * pi / per * yrs)

# Add them together so we have signal and noise
x1 <- x + wav

# Here is the redfit spec
redf.x <- redfit(x1, nsim = 500)

op <- par(no.readonly = TRUE) # Save to reset on exit
par(tcl = 0.5, mar = rep(2.2, 4), mgp = c(1.1, 0.1, 0))

plot(redf.x[["freq"]], redf.x[["gxxc"]],
     ylim = range(redf.x[["ci99"]], redf.x[["gxxc"]]),
     type = "n", ylab = "Spectrum (dB)", xlab = "Frequency (1/yr)",
     axes = FALSE)
grid()
lines(redf.x[["freq"]], redf.x[["gxxc"]], col = "#1B9E77")
lines(redf.x[["freq"]], redf.x[["ci99"]], col = "#D95F02")
lines(redf.x[["freq"]], redf.x[["ci95"]], col = "#7570B3")
lines(redf.x[["freq"]], redf.x[["ci90"]], col = "#E7298A")
freqs <- pretty(redf.x[["freq"]])
pers <- round(1 / freqs, 2)
axis(1, at = freqs, labels = TRUE)
axis(3, at = freqs, labels = pers)
mtext(text = "Period (yr)", side = 3, line = 1.1)
axis(2); axis(4)
legend("topright", c("x", "CI99", "CI95", "CI90"), lwd = 2,
       col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"),
       bg = "white")
box()
par(op)


### Superimposed epoch analysis: for testing whether respons to event years was significnatly different from background
library(graphics)
data(cana157)
event.years <- c(1631, 1742, 1845)
cana157.sea <- sea(cana157, event.years)
foo <- cana157.sea$se.unscaled
names(foo) <- cana157.sea$lag
barplot(foo, col = ifelse(cana157.sea$p < 0.05, "grey30", "grey75"), 
        ylab = "RWI", xlab = "Superposed Epoch")