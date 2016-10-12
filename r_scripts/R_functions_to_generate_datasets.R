# Function to generate datasets

is.installed = function(pkg){ 
  # Function to check whether package is installed or needs to be installed
  #   from: http://stackoverflow.com/questions/9341635/check-for-installed-packages-before-running-install-packages  
  # 
  # Arg: 
  #   pkg: package name
  # Return: NULL. Ps: named package is installed case it isn't in the system
  is.element(pkg, installed.packages()[,1])
} 

props = function(ncol, nrow, var.names=NULL){
  # Function to generate random proportions whose rowSums = 1
  # 	from: https://www.r-bloggers.com/function-to-generate-a-random-data-set/
  #
  # Args:
  #   ncol: number of columns in the dataset
  #	  nrow: number of rows in the dataset
  #   var.names: names of the variables (cols) in the dataset
  #    	         if set to NULL, it gives the default col names
  #    	         if user set variable names - e.g. c('blue', 'red', 'green') this
  #    	           should at the same size of ncols
  #
  # Returns:
  #   A data frame comprised of ncol, nrow, and var.names

# Error handling
    if (ncol < 2) stop("ncol must be greater than 1")

    p = function(n){
        y = 0
        z = sapply(seq_len(n-1), 
        			
        			function(i) {
                			x = sample(seq(0, 1-y, by=0.01), 1)
                			y = y + x
                			
                			return(x)
               				}
               		)

        w = c(z , 1-sum(z))
        
        return(w)
    }

    df = data.frame(t(replicate(nrow, p(n=ncol))))
    
    if (!is.null(var.names)) {
    	colnames(DF) <- var.names
    	}

    return(df)
}


NAinsert <- function(df, prop=0.1){
  # Function to randomly insert a certain proportion of NAs into a dataframe
  # 	from: https://www.r-bloggers.com/function-to-generate-a-random-data-set/
  #
  # Args:
  #  df: dataframe of any kind of data
  #  prop: proportion of NA values to be inserted into the data frame     
  # Returns:
  #   same input data frame but if random 'NaNs' added to it
  
  # Error handling
	if (is.data.frame(df) == FALSE) stop("df must be a data frame") 

    n = nrow(df)
    m = ncol(df)
    num.to.na = ceiling(prop*n*m)
    
    id = sample(0:(m*n-1), num.to.na, replace = FALSE)
    rows = id %/% m+1
    cols = id %% m+1
    sapply(seq(num.to.na), 

    	function(x){
          df[rows[x], cols[x]] <<- NA
        }
    )
    return(df)
}


DFmaker <- function(n=10, type=wide, digits=2, na.rate=0.0){
  # Function to randomly generate an n-lenght data set with predefined variables
  # 	from: https://www.r-bloggers.com/function-to-generate-a-random-data-set/
  # 
  # Args:
  #   n: number of id rows
  #   type: argument (default “wide” or “long”)
  #         wide = number of rows is equal the argument 'n'
  #         long = number of rows is three times bigger than the argument 'n'
  #   digits: number of decimal digits
  #   na.rate: (a decimal value between 0 and 1; default is 0) that randomly inserts missing data 
  #             great for teaching demos and testing corner cases
  #             PS: To work, it requires loading the function NAinsert
  # Return:
  #   data frame with the following col: id, group, hs.grad, race, gender, age, m.status,
  #                                      political, n.kids, income, score, time1, time2, 
  #                                      time3 (time are cumulative sums)
  
    rownamer = function(dataframe){
                  x = as.data.frame(dataframe)
                  rownames(x) <- NULL
                
                return(x)
                }

    dfround = function(dataframe, digits = 0){
                df = dataframe
                df[,sapply(df, is.numeric)] <- round(df[,sapply(df, is.numeric)], digits) 
      
               return(df)
               }

    TYPE = as.character(substitute(type))
    time1 = sample(1:100, n, replace = TRUE) + abs(rnorm(n))
    
    DF <- data.frame(id = paste0("ID.", 1:n), 
                     group = sample(c("control", "treat"), n, replace = TRUE),
                     hs.grad = sample(c("yes", "no"), n, replace = TRUE), 
                     race = sample(c("black", "white", "asian"), n, replace = TRUE, 
                                  prob=c(.25, .5, .25)), 
                     gender = sample(c("male", "female"), n, replace = TRUE), 
                     age = sample(18:40, n, replace = TRUE),
                     m.status = sample(c("never", "married", "divorced", "widowed"), 
                                       n, replace = TRUE, prob=c(.25, .4, .3, .05)), 
                     political = sample(c("democrat", "republican", "independent", 
                                          "other"), n, replace= TRUE, 
                                           prob=c(.35, .35, .20, .1)),
                     n.kids = rpois(n, 1.5), 
                     income = sample(c(seq(0, 30000, by=1000), seq(0, 150000, by=1000)), 
                                     n, replace=TRUE),
                     score = rnorm(n), time1, 
                     time2 = c(time1 + 2 * abs(rnorm(n))), 
                     time3 = c(time1 + (4 * abs(rnorm(n)))))

    if (na.rate!=0) {  
        DF <- NAinsert(DF, na.rate)
        }
    
    DF <- switch(TYPE, 
                wide    = DF, 
                long    = {DF <- reshape(DF, direction = "long", idvar = "id",
                varying = c("time1","time2", "time3"),
                v.names = c("value"),
                timevar = "time", times = c("time1", "time2", "time3"))
            rownamer(DF)}, 
            stop("Invalid Data \"type\""))

    return(dfround(DF, digits=digits))
}


regCoordinates = function(N = 50, sample = 100){
  # Function from Package ‘geosphere’ 
  #      docs at: https://cran.r-project.org/web/packages/geosphere/geosphere.pdf
  #       change: add function 'sample' to get a sample of GPS coords
  #             : add paste0 to include id names
  #             : add data.frame to create a data frame
  # Args:
  #   N: Number of ’parts’ in which the earth is subdived
  #   sample: number of samples to be drawn 
  # Return:
  #   data frame with the follwing cols: ID, log, lat, group  
  
  N = round(N)
  sample = round(sample)
  
  # Error handling 
  if (N < 1 | sample < 3) {
    stop("N should be >= 1 and n should be >= 3")
  }
  
  # Random GPS coords
  beta = 0.5 * pi/N
  A = 2 * sin(beta/2)
  points = rbind(c(0, 0, 1), c(0, 0, -1))
  R = sin(1:N * beta)
  Z = cos(1:N * beta)
  M = round(R * 2 * pi/A)
  
  for (i in 1:N) {
    j = 0:(M[i] - 1)
    Alpha = j/M[i] * 2 * pi
    X = cos(Alpha) * R[i]
    Y = sin(Alpha) * R[i]
    points = rbind(points, cbind(X, Y, Z[i]))
    
    if (i != N) {
      points = rbind(points, cbind(X, Y, -Z[i]))
    }
  }
  
  r = sqrt(points[, 1]^2 + points[, 2]^2 + points[, 3]^2)
  theta = acos(points[, 3]/r)
  phi = atan2(points[, 2], points[, 1])
  lat = theta * 180/pi - 90
  lon = phi * 180/pi
  
  lat_sample = sample(lat, sample, replace = F)
  log_sample = sample(lon, sample, replace = F)
  
  # Data frame
  DF = data.frame(id  = paste0("ID.", 1:sample),
                  log = log_sample,
                  lat = lat_sample,
                  groups = sample(c("A", "B", "C"), sample, replace = TRUE)
                  )
                  
  return(DF)
}