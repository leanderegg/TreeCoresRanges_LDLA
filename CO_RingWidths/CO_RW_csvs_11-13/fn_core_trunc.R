#######################################
#        core_trunc()
#######################################

# Use: takes a data set produced by core_average()
  # and turns cores into columns and years into rows,
  # with Elev-Tree is colnames
  # also truncates to desired length of time

# Inupts: workingdata= data frame produced by core_average
        # first = year to start time series, as character
        # last = year to end time series, as character
        # detmethod = method (default "Spline") with which to detrend
        # detrend = logical, whether to return detrended time series
            #       or undetrended ring widths, default = TRUE

#______________________________________
core_trunc <- function (workingdata, first=start, last=end, detmethod="Spline", detrend=TRUE, drop=3, yrcored=2013){
require(dplR)

# Use: takes a data set produced by core_average()
# and turns cores into columns and years into rows,
# with Elev-Tree is colnames
# also truncates to desired length of time
 if(is.na(suppressWarnings(as.numeric(colnames(workingdata)[ncol(workingdata)])))){
   colnames(workingdata)[6:ncol(workingdata)]<- seq(yrcored,yrcored-(ncol(workingdata)-6), by=-1)
 }
  # grabbing desired rows + drop from dataframe, unless there aren't enough rows to begin with
  First <- as.character(max((as.numeric(first) - drop), as.numeric(colnames(workingdata)[ncol(workingdata)])))
  if (as.numeric(First)>as.numeric(first)) print ("Time series shorter than desired")  
# establish vector of column names
  yrs <- seq(from=which(colnames(workingdata)==last),to=which(colnames(workingdata)==First), by=1)
  
  # subsetting working data to just desired range
  treedata2011 <- workingdata[,c(1:4,yrs)]
  
  treedata2011$ID <- paste(treedata2011$Transect,treedata2011$Species,treedata2011$Elev, treedata2011$Tree, sep = "-")
  idcol <- as.numeric(ncol(treedata2011))
  treedata<- treedata2011[,c(1:4,idcol,5:(idcol-1))]  # good to go!
  #________________________________________________________________
  
  
  
  
  #------------------------------------------------------------
  ## Creating a detrended rw data frame (lifted from Ailene's xdating code)
  lastyr <- which(colnames(treedata)==last)
  spstdat.rwl<-t(treedata[,lastyr:ncol(treedata)]) ### 
  rownames(spstdat.rwl)<-colnames(treedata[lastyr:ncol(treedata)])
  colnames(spstdat.rwl)<-treedata$ID
  spstdat.rwl.df<-as.data.frame(spstdat.rwl)
  
if (drop > 0){
  for (i in 1:ncol(spstdat.rwl.df)){
    tmp <- spstdat.rwl.df[,i]
    if (is.na(spstdat.rwl.df[length(tmp),i])){
      l <- min(which(is.na(tmp)))
      drops <- seq(from=l-drop, to=l-1)
      spstdat.rwl.df[drops,i] <- NA
    }
  }
  spstdat.drop <- spstdat.rwl.df[1:which(rownames(spstdat.rwl.df)==first),]
}
else {spstdat.drop <- spstdat.rwl.df}
  #__________ CHANGE FOR FINE TUNING DETRENDING METHOD________
  rw.det <- detrend(rwl=spstdat.drop, method=detmethod) # final object if detrending
  rw.nodet <- spstdat.drop # final object if not detrending
  #-----------------------------------------------------------
  #==========================================================
  
  if (detrend==TRUE)
    return(rw.det)
  else
    return(rw.nodet)
}