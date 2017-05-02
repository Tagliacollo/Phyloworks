clean.coords <- function(DF, ext="p")
{
	rst  <- raster(dem,xmn=-180, xmx=180, ymn=-60, ymx=90)
    data(world)
    
    ask = ""
    while (ask != "NO") 
    {
    ssp    <- DF$Species
    x      <- coord2numeric(DF$x)
    y      <- coord2numeric(DF$y)
    id     <- DF$ID
    source <- DF$source
    
    cords  <- data.frame(DF$x, DF$y)
    xy     <- SpatialPoints(cords)
    ex     <- getextent(x, y, ext)
    xlm    <- ex$xlm
    ylm    <- ex$ylm

    sp::plot(world, border = "gray", 
    	     xlim=xlm, ylim=ylm, 
    	     main=as.character(ssp[1]), axes = T)
    vals <-  extract(rst, xy)
    points(x, y, col=as.factor(source))
    legend('topright', levels(DF$source), 
    lty=1, col=1:3, bty='n', cex=.75)

    gps_point <- identify(x, y, n = 1, plot = FALSE)
    ssp_point    <- DF$Species[gps_point]
    id_point     <- DF$ID[gps_point]
    x_point      <- DF$x[gps_point]
    y_point      <- DF$y[gps_point]
	source_point <- DF$source[gps_point] 

  	text(x_point, y_point, labels = source_point, pos = 4, adj = 1)        
    mtext(paste("x = ", x_point, "; y = ", y_point, sep = ""), 
    	  			side = 4, line = 1, cex = 0.9)

    print("What would you like to do? SKIP, DELETE, or replace species name (ex. Gymnotus tigre)")
    TODO <- readline(prompt="your answer: " )

    if (TODO == "DELETE") { DF <- subset(DF, ID!=id_point) } 
    if (TODO != "SKIP") {DF <- within(DF, Species[Species==ssp_point & ID==id_point] <- TODO) }

    print("Do you still have to make changes?")
    ask <- readline(prompt="YES/NO: " )
    }

  return(DF)
}










