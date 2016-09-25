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

# Function to compare phy_nwklo vs dataset and clean up missing data
clean_up_data = function(phy_nwk, dat_csv) {
  
  phy_nwkOnly <- setdiff(phy_nwk$tip.label, rownames(dat_csv))
  datOnly <- setdiff(rownames(dat_csv), phy_nwk$tip.label)
  
  if (length(phy_nwkOnly) > 0) {
    phy_nwk <- drop.tip(phy_nwk, phy_nwkOnly)
    message("The following taxa were deleted from the tree: "); print(phy_nwkOnly)
  }
  
  if (length(datOnly) > 0) {	
    dat_csv <- dat_csv[-match(datOnly, rownames(dat_csv)),]
    message("The following taxa were deleted from the dataset: "); print(datOnly) 
  }
  
  dat_csv <- dat_csv[phy_nwk$tip.label,]
  if (rownames(dat_csv) == phy_nwk$tip.label) {
    return (list(phy_nwk, dat_csv))
  }
  #TODO else to catch error if rownames(dat_csv) != phy_nwk$tip.label
}

# Function to plot biogeo data on trees
plot_biogeo = function (phy_nwk, dat_csv) {
 
    # TODO: add the path to folders 
    pdf('plot.pdf', width=8, height=10, paper = 'a4')
    plot(phy_nwk, cex=0.5, label.offset = (ncol(dat_csv) + 1), no.margin = TRUE)
    eixo_x = max(nodeHeights(phy_nwk))
    eixo_y = length(phy_nwk$tip.label)
    
    for (i in 1:ncol(dat_csv)) {
      dat_col = dat_csv[i]
      label = dat_col[match(phy_nwk$tip.label, row.names(dat_col)), ]
      
      points(rep(eixo_x + i, eixo_y), 1:(eixo_y), pch=22, bg = label, cex=1)
      text((eixo_x + i), (eixo_y + 0.8), colnames(dat_col), cex=0.4, srt=90, adj=0.0)
    }
    
    dev.off()
    return (list(phy_nwk, dat_csv))
}
