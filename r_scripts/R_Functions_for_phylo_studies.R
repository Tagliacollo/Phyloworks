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

#Function to get best model of discrete character evolution
eval_fitDiscrete_models = function (phy_nwk, dat_csv, numcol = 1){
  
  dat_vector = as.factor(dat_csv[, numcol]) 
  names(dat_vector) = row.names(dat_csv)
  
  model_choices = c('ER', 'SYM', 'ARD')
    
  fit_model = list()
  for (i in 1:3) {
    fit = fitDiscrete(phy_nwk, dat_vector, model = model_choices[i])
    fit_model[[model_choices[i]]] = fit
  }
  
  #csv table
  lnL      = c(fit_model$ER$opt$lnL, fit_model$SYM$opt$lnL, fit_model$ARD$opt$lnL)
  aicc     = c(fit_model$ER$opt$aicc, fit_model$SYM$opt$aicc, fit_model$ARD$opt$aicc)
  del_aicc = c(aicw(aicc)$delta)
  wts_aicc = c(aicw(aicc)$w)
  
  output_csv = matrix( , 3, 4, dimnames=list(c("ER", "SYM", "ARD"),
                                             c("lnL", "AICc", "delta AICc", "AICc weights")))
  output_csv[,1] = lnL  
  output_csv[,2] = aicc
  output_csv[,3] = del_aicc
  output_csv[,4] = wts_aicc
  
  write.csv(output_csv, file = 'model_selection.csv')
  
  best_fitted_model = fit_model[which.min(aicc)]
  
  return (list(best_fitted_model, dat_vector))
}

# Fuction to make simmaps and plot 1st simmap result
plot_simmmap_pdf = function (phy_nwk, dat_vector, model, nsim = 100){
  mtrees = make.simmap(phy_nwk, dat_vector, model, nsim)
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