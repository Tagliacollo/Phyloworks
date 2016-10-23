l1ou.models = function(phy, dat, criterion = 'AICc', name = 'l1ou_models.pdf'){
  
  dat.mtx = as.matrix(dat)
  dat.processed = adjust_data(phy, dat.mtx)
  
  # regimes
  shift.model = estimate_shift_configuration(dat.processed$tree, dat.processed$Y)
  conv.regimes = estimate_convergent_regimes(shift.model, criterion = criterion)
  
  return(shift.model = shift.model, conv.regimes = conv.regimes)
  
  models.l1ou = list(shift.model, conv.regimes)
  
  # save pdf
  pdf(name, width=8, height=11)
  
  # plot regimes
  for (model in models.l1ou) {
    plot(model, show.tip.label = FALSE,  asterisk = F, edge.ann.cex = 0.6)
    
  }
  
  dev.off()
  
  cmdstr = paste("open ", name, sep="")
  system(cmdstr)
  
}

l1ou.bootstrap = function(shift.model, conv.regimes, n.bootstrap, name = 'l1ou_boots.pdf', 
                          multicore = FALSE, nCores = 2, quietly = TRUE ) {
  
  shift.model.boots = l1ou_bootstrap_support(shift.model, nItrs = n.bootstrap, 
                                             multicore, nCores, quietly)
  conv.model.boots = l1ou_bootstrap_support(conv.regimes, nItrs = n.bootstrap, 
                                            multicore, nCores, quietly)
  
  return(shift.model.boots = shift.model.boots, conv.regimes.boots = conv.regimes.boots)
  
  res.l1ou.boots = list(shift.model.boots, conv.model.boots)
  models.l1ou = list(shift.model, conv.regimes)
  
  # save pdf
  pdf(name, width=8, height=11)
  
  # plot regimes
  for (boots in res.l1ou.boots) {
    nodes <- round(boots$detection.rate*100, digits=1)
    nodes <- ifelse(nodes > 50, paste0(nodes), NA)
    
    i = 1
    plot(models.l1ou[[i]], edge.label = nodes, edge.ann.cex = 0.5, edge.label.ann = TRUE, 
         show.tip.label = FALSE,  edge.shift.ann = F, asterisk = F, 
         edge.label.adj = c(1.0, 1.0))
    i = i + 1
  }
  
  dev.off()
  
  cmdstr = paste("open ", name, sep="")
  system(cmdstr)
}
