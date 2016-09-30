library(phytools)

pw_plot_Phylomorphospace = function(phy, dat.csv, anc.estimates = FALSE, 
                                             pdf.name = NULL, ...) {
  # Plot a phylomorphoscape with information of ancestral area estimates
  #
  # Args:
  #  phy:  Phylogeny stored as the class "phylo"
  #  dat.csv:  CSV file with species names as rownames.
  #            Be sure your original CSV has four columns being the 1º comprised of species 
  #              name, 2º of a discrete character and 3º e 4º of lat and lon values (in decimals)
  #            Be sure phy and dat.csv have same length and are in the same order.
  #              Use PwCleanUp to prepare the input data
  #  anc.estimates: estimates of ancestral states obtained from Ace (package ape)
  #                 if anc.estimates is set to TRUE, it plots estimates of ancestral states 
  #                 at each node of the phylo
  #  pdf.name: name of pdf plot saved in the working directory 
  # Returns:
  #  NULL. PDF plot is saved in the working directory
  
  # Error handling
  if (rownames(dat.csv) != phy$tip.label){
    stop('rownames of dat.csv and tip.labels of phy do not match each other. Use PwCleanUp!')
  }
  
  if (ncol(dat.csv) != 3) {
    stop('dat.csv must have 3 data columns')
  }
  
  if (is.factor(dat.csv[,1]) == FALSE) {
    stop('The columns must contains a discrete (str: factor) and two continuous characters, respectively')
  }
  
  if (is.null(pdf.name)) {
    stop('provide a pdf.name to your plot - pdf.name = NULL')
  }
#  if (anc.estimates  ) {
#    stop('provide a ** ')
  
  # Vectors
  vec.discrete   = as.factor(dat.csv[,1]);  names(vec.discrete)   = row.names(dat.csv) 
  vec.continuous = dat.csv[,2:3];           # colnames(vec.continuous) = colnames(dat.csv) 
  
  # Plot parameters
  cols = setNames(rainbow(nlevels(vec.discrete), alpha = 1.0)
                  [1:length(unique(vec.discrete))], sort(unique(vec.discrete)))
  
  col.tips = as.character(vec.discrete);   names(col.tips) <- phy$tip.label
  
  # Plot - save 'pdf.name' in the working directory
  name = paste(pdf.name, '.pdf', sep = '')
  pdf(name, width=7, height=7)
  
  phylomorphospace(phy, vec.continuous, ftype="off", node.size=rep(0,2))
  tiplabels(pie=to.matrix(col.tips, levels(vec.discrete)), piecol=cols, cex= 0.5)
  
  legend("topleft", legend=names(cols), col=cols, fill = cols, bty = 'n')
         
  if (anc.estimates == TRUE) {
    lik.anc = ace(col.tips, phy, type="discrete", model="ER")
    nodelabels(pie=lik.anc$lik.anc, piecol=cols, cex=0.6)
  }
  dev.off()
  
  cmdstr = paste("open ", name, sep="")
  system(cmdstr)
  }


pw_plot_BiogeoAreas = function(phy, dat.csv, pdf.name = NULL, ...) {
  # Plot species area distributions on a phylo
  #
  # Args:
  #  phy:  Phylogeny stored as the class "phylo"
  #  dat.csv:  CSV file with species names as rownames and areas as colnames.
  #            Be sure phy and dat.csv have same length and are in the same order.
  #            Use PwCleanUp to prepare the input data
  #  pdf.name: name of pdf plot saved in the working directory 
  # Returns:
  #  NULL. PDF plot is saved in the working directory
  
  # Error handling
  if (rownames(dat.csv) != phy$tip.label){
    stop('rownames of dat.csv and tip.labels of phy do not match. Use PwCleanUp!')
  }
  
  if (is.null(pdf.name)){
    stop('provide a pdf.name to your plot - pdf.name = NULL')
  }
  
  # Plot parameters
  label.offset = (max(nodeHeights(phy)) / 50) * ncol(dat.csv) 
  eixo.x = max(nodeHeights(phy))
  eixo.y = length(phy$tip.label)
  pont.txt = label.offset / ncol(dat.csv)  
  
  # Plot - save 'pdf.name' in the working directory
  name = paste(pdf.name, '.pdf', sep = '')
  pdf(name, width=8, height=10, paper = 'a4')
  
  plot(phy, cex=0.5, label.offset = label.offset, no.margin = TRUE)
  add.scale.bar()
  
  j = 0 # j is the counter for plotting points and texts in the x.axis
  for (i in 1:ncol(dat.csv)) {
    dat.col = dat.csv[i]
    label = dat.col[match(phy$tip.label, row.names(dat.col)), ]
    
    points(rep(eixo.x + j, eixo.y), 1:(eixo.y), pch=22, bg = label, cex=1.7)
    text((eixo.x + j), (eixo.y + 0.8), colnames(dat.col), cex=0.4, srt=90, adj=0.0)
    j = j + pont.txt
    
  }
  dev.off()
  
  # open plot in a PDF reader
  cmdstr = paste("open ", name, sep="")
  system(cmdstr)
}