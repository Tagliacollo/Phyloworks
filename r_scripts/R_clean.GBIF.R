clean.GBIF <- function(DF, vars=c("bio1", "bio12")) 
{

  require(biogeo)

  DF <- addmainfields(DF, species="Species")
  checkdatastr(DF)
  
  print(" ... cleaning up GBIF coordinates, please wait")

  dem  <- raster(dem,xmn=-180, xmx=180, ymn=-60, ymx=90)
  data(world)

  # check this post if you get the following msg: NOTE: rgdal::checkCRSArgs: no proj_defs.dat in PROJ.4 shared files
  # https://gis.stackexchange.com/questions/224467/raster-error-note-rgdalcheckcrsargs-no-proj-defs-dat-in-proj-4-shared-file/224700

  cleanup  <- errorcheck(world, dem, DF, countries="country",
                         countryfield="NAME", vars=vars, res=10, 
                         elevc="elevation", diff=50)

# DF <- cleanup[cleanup$error!=1 & cleanup$dups!=1, ]
  DF <- cleanup[cleanup$dups!=1 & cleanup$wrongEnv!=1, ]
 
  return(DF[rowSums(is.na(DF)) != ncol(DF),])
}