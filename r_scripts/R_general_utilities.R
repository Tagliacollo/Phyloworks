is_installed = function(pkg)
  { 
  # Function to check whether package is installed or needs to be installed
  #   from: http://stackoverflow.com/questions/9341635/check-for-installed-packages-before-running-install-packages  
  # 
  # Arg: 
  #   pkg: package name
  # Return: NULL. Ps: named package is installed case it isn't in the system
  is.element(pkg, installed.packages()[,1])
  } 

source_https = function(url, ...) 
  {
  # Function to source R scripts hosted on github repository
  #   from: https://tonybreyal.wordpress.com/2011/11/24/source_https-sourcing-an-r-script-from-github/
  
  # Arg: 
  #   url: https github path (raw path)
  # Return: all functions in the script
  
  require(RCurl)
  
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) 
    {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
    })
  }


set_dir_path = function(dataset.path, parent.name = 'processed_data') 
  {
  # Function to set working directory in R
  # Arg: 
  #   dataset.path: full dataset path to a csv file
  #   parent.name: name of the parent folder to be added to the base path
  #
  # Return: set new working directory

  base   = dirname(dirname(dataset.path))
  parent = paste(base, '/', parent.name, sep = "") 
  child  = strsplit(basename(dataset.path), "\\.")[[1]][1]
  
  path = dir.create(file.path(parent, child))
  cat('dir path set to',  file.path(parent, child))
  
  return(setwd(file.path(parent, child)))
  }