ggmap.GGMAP <- function(DF) 
{
	require(ggmap)

	DF <- data.frame(DF)
	sp_name <- as.vector(DF$Species[1])

	lat <- c(min(DF$y), max(DF$y))
	lon <- c(min(DF$x), max(DF$x))

	myMap  <- get_map(location=c(mean(lon), mean(lat)), 
		              zoom = 4, maptype ="satellite", source = "google", 
		              crop=FALSE)
	ggmap(myMap) +
	geom_point(aes(x=DF$x, y=DF$y), data=DF,
     	       alpha=0.6, color="darkred", size=1) +
			   ggtitle(sp_name)
}