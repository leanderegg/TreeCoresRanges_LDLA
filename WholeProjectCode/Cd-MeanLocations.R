###############################################
#          Transect Locations
#         And plotting code for catching Lat/Long screw ups
################################################
# takes tree location data from Locations.csv and produces mean locations
# code plotting climateWNA normals has been outsourced to Cd-Fig-FlimateNormals.R in the 'climate vis' folder
# as of 11/1/15

# now instead contains code for plotting location data to verify it.

require(RColorBrewer)
require(scatterplot3d)
require(rgl)
require(Rcmdr)
require(reshape)
require(dplyr)

locs <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WholeProjectCode/Locations.csv", header=T)
# states <- unique(locs$State)
# 
# mean.locs <- c()
# for (i in states){
#   tmp.state <- subset(x=locs,subset=State==i)
#   for (j in unique(tmp.state$Species)){
#     tmp.spp <- subset(tmp.state, subset= Species==j)
#     for (k in unique(tmp.spp$Band)){
#       tmp.band <- subset(tmp.spp, subset= Band==k)
#       mLat <- mean(tmp.band$Lat)
#       mLon <- mean(tmp.band$Lon)
#       mElev <- mean(tmp.band$Elevation)
#       mean.locs<- rbind(mean.locs, c(i,j,k, mLat, mLon, mElev))
#     }
#     
#   }
# }
# colnames(mean.locs) <- c("State", "Species", "Band", "mLat", "mLon", "mElev")
# 
# m.lat <- tapply( locs$Lat,INDEX=list(locs$Band, locs$Species, locs$State), FUN=mean, simplify=T)
# m.lat[,,1]
# dim(m.lat)




meanlocs <- locs %>% group_by (State, Species, Band) %>% summarise (mLat = mean(Lat), mLon=mean(Lon), mElev = mean(Elevation))
# write.csv(meanlocs, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WholeProjectCode//MeanLocations_7_31_15.csv")

wnalocs <- cbind(as.character(meanlocs$State), paste(meanlocs$Species,meanlocs$Band, sep="-"), meanlocs$mLat, meanlocs$mLon, meanlocs$mElev)
colnames(wnalocs) <- c("ID1", "ID2", "lat", "lon","el")
write.csv(meanlocs, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WholeProjectCode//forclimateWNA_MeanLocations_8_5_15.csv")



############################## *** MAPPING to check location data **** ##################
library(RgoogleMaps)

temp <- locs[which(locs$Species=="PIPO"),]
plot(temp$Lat~I(-1*temp$Lon), col=temp$Band)
temp <- locs[which(locs$Species=="POTR"),]
plot(temp$Lat~I(-1*temp$Lon), col=temp$Band)
temp <- locs[which(locs$Species=="ABLA"&locs$State=="CO"),]
plot(temp$Lat~I(-1*temp$Lon), col=temp$Band)

temp <- locs[which(locs$Species=="TSHE"&locs$State=="MT"),]
plot(temp$Lat~I(-1*temp$Lon), col=temp$Band)
temp <- locs[which(locs$Species=="PSME"&locs$State=="MT"),]
plot(temp$Lat~I(-1*temp$Lon), col=temp$Band)
temp <- locs[which(locs$Species=="ABLA"&locs$State=="MT"),]
plot(temp$Lat~I(-1*temp$Lon), col=temp$Band)

temp <- locs[which(locs$Species=="TSHE"&locs$State=="WA"),]
plot(temp$Lat~I(-1*temp$Lon), col=temp$Band)
temp <- locs[which(locs$Species=="PSME"&locs$State=="WA"),]
plot(temp$Lat~I(-1*temp$Lon), col=temp$Band)
temp <- locs[which(locs$Species=="ABLA"&locs$State=="WA"),]
plot(temp$Lat~I(-1*temp$Lon), col=temp$Band)

locs <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WholeProjectCode/Locations.csv", header=T)


## getting the appropriate google map
siteTerrain <- GetMap(center=c(lat=mean(temp$Lat), lon=-1*mean(temp$Lon)), #37.47,lon=-108.24),
                      destfile="/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/Maps/tmpWAmap.png",
                      maptype="hybrid",zoom=15)

## adding missing trees as pts
tmp <- PlotOnStaticMap(siteTerrain, lat = temp$Lat,
                       lon = -1*temp$Lon, 
                       destfile = "RgoogleMap3.png", cex=1.5,pch=20,
                       col="green",FUN=points,add=F)
## adding tree names of missing trees
TextOnStaticMap(siteTerrain, lat = temp$Lat,
                lon = I(-1*temp$Lon), 
                cex=1,labels=temp$Pair,col="red", add=T);
## adding 2013 pair names to narrow it down.
TextOnStaticMap(siteTerrain, lat = trees13$Lat,
                lon = I(-1*trees13$Long), 
                cex=1,labels=trees13$Pair,col="red", add=T);

