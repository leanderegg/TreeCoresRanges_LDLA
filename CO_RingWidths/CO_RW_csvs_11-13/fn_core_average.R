####### Functions for taking WinDendro output and cleaning it up into just ring width dataframes


core_average <- function(data, yrcored = 2013)
{ 
  ## This function takes a raw csv from crossdating with dplR and WinDendro,
  # pulls out most recent version of a chronology, removes all the Windendro output columns
  # and then averages the cores from each tree, pulling out NAs.
data$Analysis.Date.Time <- as.POSIXct(strptime(data$Analysis.Date.Time, format = "%m/%d/%Y %H:%M"))
un <- unique(data$TreeID)
dataclean <- data[1:length(un),]
for(i in 1:length(un))
{
data_subset <- data[which(data$TreeID==un[i]),]
dataclean[i,] <- data_subset[which(data_subset$Analysis.Date.Time  == max(data_subset$Analysis.Date.Time)),]
print(i)
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
#coredata[,"X"] <- NULL # removing 2013 (year core was collected)

trees <- unique(coredata$Tree)
Transect <- rep(coredata$Transect[1], times=length(trees))
Species <- rep(coredata$Species[1], times=length(trees))
Elev <- rep(coredata$Elev[1], times=length(trees))
averages <- aggregate(x=coredata[,10:ncol(coredata)],by=list(coredata$Tree),FUN="mean", na.rm=T)
nyears <- ncol(averages)-2
colnames(averages)<- c("Tree", seq(from=yrcored, to= yrcored-nyears, by=-1))
treedata <- cbind(Transect,Species,Elev,averages)

return(treedata)
}



#############################
#     Function to keep each core unaveraged
#############################
# essentially the same as core_average(), but doesn't average cores per tree
core_clean <- function(data, yrcored = 2013)
{
  # essentially the same as core_average(), but doesn't average cores per tree
  data$Analysis.Date.Time <- as.POSIXct(strptime(data$Analysis.Date.Time, format = "%m/%d/%Y %H:%M"))
  un <- unique(data$TreeID)
  dataclean <- data[1:length(un),]
  for(i in 1:length(un))
  {
    data_subset <- data[which(data$TreeID==un[i]),]
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
  #coredata[,"X"] <- NULL # removing 2013 (year core was collected)
  #  return(coredata)
  trees <- unique(coredata$TreeID)
  Transect <- rep(coredata$Transect[1], times=length(trees))
  Species <- rep(coredata$Species[1], times=length(trees))
  Elev <- rep(coredata$Elev[1], times=length(trees))
  # replaced $Tree with $TreeID. Theoretically this keeps all cores.
  averages <- aggregate(x=coredata[,10:ncol(coredata)],by=list(coredata$TreeID),FUN="mean", na.rm=T)
  nyears <- ncol(averages)-2
  colnames(averages)<- c("CoreID", seq(from=yrcored, to= yrcored-nyears, by=-1))
  treedata <- cbind(Transect,Species,Elev,averages)
  
  return(treedata)
}

