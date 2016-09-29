library(ape)
library(phytools)
library(geiger)

# Function to rescale & ladderize phylogenies 
rescale_phy = function(phy, scale = 1, ladderize = TRUE) {
  #From: http://blog.phytools.org/2012/02/quicker-way-to-rescale-total-length-of.html
  phy$edge.length = phy$edge.length/max(nodeHeights(phy)[,2])*scale
  
  if (ladderize == TRUE) {
    phy_ladderize = ladderize(phy)
    return (phy_ladderize)
  }
  else {
    return (phy)
  } 
}

# Function to compare phylo vs dataset and clean up missing data
pw_Cleanup_data = function(phy, dat.csv) {
  
  dat.csv <- dat.csv[phy$tip.label,]
  
  if (rownames(dat.csv) == phy$tip.label) {
    } else {
      
      #add info about read.csv, rownames = TRUE
      phyOnly <- setdiff(phy$tip.label, rownames(dat.csv))
      datOnly <- setdiff(rownames(dat.csv), phy$tip.label)
      
      if (length(phyOnly) > 0) {
        phy <- drop.tip(phy, phyOnly)
        message("The following taxa were deleted from the tree: "); print(phyOnly)
        }
      
      if (length(datOnly) > 0) {
        dat.csv <- dat.csv[-match(datOnly, rownames(dat.csv)),]
        message("The following taxa were deleted from the dataset: "); print(datOnly)
        }
      }
  message('You are ready to proceed! see phy: [[1]], csv: [[2]]')
  return (list(phy, dat.csv))
  }

#Function to get best model of discrete character evolution
pw_fitDiscreteModels = function (phy, dat.csv, numcol = 1){
  
  dat.vector = as.factor(dat.csv[, numcol]) 
  names(dat.vector) = row.names(dat.csv)
  
  model.choices = c('ER', 'SYM', 'ARD')
    
  fit.model = list()
  for (i in 1:3) {
    fit = fitDiscrete(phy, dat.vector, model = model.choices[i])
    fit.model[[model.choices[i]]] = fit
  }
  
  #csv table
  lnL      = c(fit.model$ER$opt$lnL, fit.model$SYM$opt$lnL, fit.model$ARD$opt$lnL)
  aicc     = c(fit.model$ER$opt$aicc, fit.model$SYM$opt$aicc, fit.model$ARD$opt$aicc)
  del.aicc = c(aicw(aicc)$delta)
  wts.aicc = c(aicw(aicc)$w)
  
  output.csv = matrix( , 3, 4, dimnames=list(c("ER", "SYM", "ARD"),
                                             c("lnL", "AICc", "delta AICc", "AICc weights")))
  output.csv[,1] = lnL  
  output.csv[,2] = aicc
  output.csv[,3] = del.aicc
  output.csv[,4] = wts.aicc
  
  write.csv(output.csv, file = 'model_selection.csv')
  
  best.fitted.model = fit.model[which.min(aicc)]
  
  return (best.fitted.model)
}

# Fuction to make simmaps and plot 1st simmap result
pw_plot_Simmmap = function (phy, dat_vector, model, nsim = 100){
  mtrees = make.simmap(phy, dat_vector, model, nsim)
  stoch_map = describe.simmap(mtrees, plot=FALSE)
  
  cols = setNames(palette()[1:length(unique(dat_vector))], sort(unique(dat_vector)))

  pdf('Simmap.pdf', width=8, height=10, paper = 'a4')
  plotSimmap(mtrees[[1]], cols, fsize=0.5, ftype="i")
  
  nodelabels(pie=stoch_map$ace, piecol=cols, cex=0.5)
  add.simmap.legend(colors=cols, prompt=FALSE, 
                    x=0.9*par()$usr[1], y=12) # TODO: improve the code for the axis
  dev.off()
  return (list(mtrees, stoch_map))
}