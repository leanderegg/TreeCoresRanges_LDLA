#############################################
#            Extracting Elev
#              from dems 
###############################################


# solutions suggested at: http://stackoverflow.com/questions/8973695/conversion-for-latitude-longitude-to-altitude-in-r

#1) Querry a website:
library(RCurl)
library(XML)

latitude <- 52.4822
longitude <- -1.8946
url <- paste(
  "http://www.earthtools.org/height",
  latitude, 
  longitude,
  sep = "/"
)
# another possible website: http://gisdata.usgs.net/xmlwebservices2/elevation_service.asmx?op=getElevation
page <- getURL(url)
ans <- xmlTreeParse(page, useInternalNodes = TRUE)
heightNode <- xpathApply(ans, "//meters")[[1]]
(height <- as.numeric(xmlValue(heightNode)))

#2) Use geonames and its associated models
require(geonames)
GNsrtm3(54.481084,-3.220625)
GNgtopo30(54.481084,-3.220625)


# suggestions from: http://stackoverflow.com/questions/21593868/extracting-elevation-from-website-for-lat-lon-points-in-australia-using-r

# 3) Using the google api
googEl <- function(locs)  {
  require(RJSONIO)
  locstring <- paste(do.call(paste, list(locs[, 2], locs[, 1], sep=',')),
                     collapse='|')
  u <- sprintf('http://maps.googleapis.com/maps/api/elevation/json?locations=%s&sensor=false',
               locstring)
  res <- fromJSON(u)
  out <- t(sapply(res[[1]], function(x) {
    c(x[['location']]['lat'], x[['location']]['lng'], 
      x['elevation'], x['resolution']) 
  }))    
  rownames(out) <- rownames(locs)
  return(out)
}

m <- matrix(c(146.9442, 146.4622, -36.0736, -36.0491), nc=2)

googEl(m)


# 4) use getData from {raster}

library(raster)
m <- data.frame(lon = c(146.9442, 146.4622), lat = c(-36.0736, -36.0491))

x <- getData('alt', country = "AUS")

cbind(m, alt = extract(x, m))
# or to interpolate rather than using nearest neighbor
cbind(m, alt = extract(x, m, method = "bilinear"))

# and raster layer is saved in object x
plot(x)
points(m)
