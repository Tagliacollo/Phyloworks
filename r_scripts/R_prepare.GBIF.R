prepare.GBIF <- function(DF)
{

  require(stringr)
  
  # check ranges of lat and log 
  DF <- subset(DF, decimalLongitude >= -180 & decimalLongitude <= 180 &
                   decimalLatitude  >= -90  & decimalLatitude  <= 90 ) 
  
  # ID re-name
  #  idnam   <- seq(1,length(DF$name), by=1)
  #  DF$key <- idnam

  # country's name
  sptstr  <- str_split_fixed(DF$country, ",", 2)
  DF$country <- sptstr[,1]

  # bioclim variable: check http://www.worldclim.org/bioclim
  ev        <- getData("worldclim", var="bio", res=10)
  ev        <- ev[[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                   12, 13, 14, 15, 16, 17, 18, 19)]]
  names(ev) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", 
                "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", 
                "bio13", "bio14", "bio15", "bio16", "bio17","bio18", "bio19")
  
  coords <- data.frame(DF$decimalLongitude, DF$decimalLatitude)
  points <- SpatialPoints(coords, proj4string = ev@crs) 
  values <- extract(ev, points) 

  DF     <- cbind.data.frame(DF, values)

  # columns names
  setnames(DF, "name", "Species")
  setnames(DF, "key", "ID")
  setnames(DF, "decimalLongitude", "x")
  setnames(DF, "decimalLatitude", "y")
  
  return(DF)
}

