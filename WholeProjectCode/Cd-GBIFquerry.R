require(rgbif)

pondikey <- name_suggest(q="Pinus ponderosa",rank="species")$key[1]
name_suggest(q="Pinus", rank="species") # there are some 100 pinus spp in gbif

occ_count(taxonKey = pondikey)#, georeferenced = TRUE)
  # 1589 georeferenced pondi records

pondidat <- occ_search(taxonKey = pondikey, continent="north_america", spatialIssues=FALSE, hasCoordinate=TRUE, return="data", limit=4000 )
gbifmap(input=pondidat,mapdatabase = 'usa')
tail(pondidat)[,1:10]



# Let's see if I can plot the Little range map shapefile.
getwd()

pondirange <- readShapeSpatial()