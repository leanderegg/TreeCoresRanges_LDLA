######## basal_area_inc() ##########
# generic function for taking output from core_averages(), plus DBH and bark depth data
# and turning it into basal area estimates
# there is also a function in dplR that does this. haven't compared the two though
# can also just return radius if calBAI=FALSE

# requires: 
# ringwidths = ring width data frame from core_average(), with $Tag
# DBHbark = tree DBH and bark info, with $Tag
# DBHcol and barkcol = column numbers where to find DBH and bark
basal_area_inc <- function(ringwidths, DBHbark, DBHcol=12, barkcol=13, calBAI=TRUE){
### Generic function for taking ring widths + DBH + barkdepth and --> BAI or radius
  # requires: 
  # ringwidths = ring width data frame from core_average(), with $Tag
  # DBHbark = tree DBH and bark info, with $Tag
  # DBHcol and barkcol = column numbers where to find DBH and bark
    
  
  #make a Tag for ringwidths that works with my DBH data frame
  #get the initial radius ()
  radinit <- ((DBHbark[,DBHcol])/2)-DBHbark[,barkcol] # in cm
  DBHbark$radinit <- radinit *10 #in mm
  # order them so that they're in the same order
  # DBHbark <- DBHbark[order(DBHbark$Tag),]
  if(length(which(DBHbark$Tag!=ringwidths$Tag))>0){
    print("order not the same in rw and tree info")
  }else{
    #ringwidths <- ringwidths[order(ringwidths$Tag),]
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
    if (calBAI==FALSE){
      # add in tree tages as row names
      rownames(radius) <- DBHbark$Tag
      # add in years as column names. but get rid of any X or X.1,X.2 at beginning from exporting csvs, and tack on a final year at the beginning
      lastyear <- paste0("X",as.numeric(strsplit(colnames(ringwidths)[ncol(ringwidths)], split="X")[[1]][2])-1)
      yearnames <- colnames(ringwidths)[grep("X", colnames(ringwidths))]
      colnames(radius) <- c(yearnames[which(nchar(yearnames)>3)],lastyear)
      return(radius)
    }
    else{
      BAI <- matrix(nrow=dim(radius)[1], ncol=dim(radius)[2])
      for(i in 1:(ncol(radius)-1)){ # this will return a matrix with 1 fewer column than radius
        BAI[,i] <- (radius[,i]^2 *pi) - (radius[,i+1]^2 *pi)
      }
      rownames(BAI) <- DBHbark$Tag
      # add in years as column names. but get rid of any X or X.1,X.2 at beginning from exporting csvs, and tack on a final year at the beginning
      lastyear <- paste0("X",as.numeric(strsplit(colnames(ringwidths)[ncol(ringwidths)], split="X")[[1]][2])-1)
      yearnames <- colnames(ringwidths)[grep("X", colnames(ringwidths))]
      colnames(BAI) <- c(yearnames[which(nchar(yearnames)>3)],lastyear)
      return(BAI)  
    }
  }}



