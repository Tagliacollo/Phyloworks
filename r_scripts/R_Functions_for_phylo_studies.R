library(ape)
library(phytools)
library(geiger)

# Function to rescale & ladderize phylogenies 
rescale_phy = function(phy_nwk, scale = 1, ladderize = TRUE) {
  #From: http://blog.phy_nwktools.org/2012/02/quicker-way-to-rescale-total-length-of.html
  phy_nwk$edge.length = phy_nwk$edge.length/max(nodeHeights(phy_nwk)[,2])*scale
  
  if (ladderize == TRUE) {
    phy_nwk_ladderize = ladderize(phy_nwk)
    return (phy_nwk_ladderize)
  }
  else {
    return (phy_nwk)
  } 
}
