library(ape)
library(l1ou)

source.https <- function(url, ...) {
  # Function to source R scripts hosted on github repository
  #   from: https://tonybreyal.wordpress.com/2011/11/24/source_https-sourcing-an-r-script-from-github/
  
  # Arg: 
  #   url: https github path (raw path)
  # Return: all functions in the script
  
  require(RCurl)
  
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}
source.https('https://raw.githubusercontent.com/Tagliacollo/Phyloworks/master/projects/KEvans/Func_conv_evol-l1ou.R')

phy = read.tree(file.choose())
dat = read.csv(file.choose(), row.names = 1)

models = l1ou.models(phy, dat, criterion = 'pBIC', name = 'kory_l1ou_models.pdf')

boots = l1ou.bootstrap(models$shift, models$convegence, 2, 
                       name = 'kory_l1ou_boots.pdf')
