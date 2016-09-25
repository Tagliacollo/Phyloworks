library(ape)
library(phytools)
library(geiger)
wd = '/Volumes/VATagliacollo/Pesquisa/Projects/Roxo_clade_HNO-Macroevolution/FFR-data'
setwd(wd)

phy = read.tree(file.choose())
dat = read.csv(file.choose(), row.names=1, header=TRUE)

clean_data = clean_up_data(phy, dat)
phy_new = clean_data[[1]]
dat_new = clean_data[[2]]


#plot discrete data
dat_plot = dat_new[-1]
plot_biogeo(phy_new, dat_plot)

#models of discrete evolution. The best model is used to estimate simmaps
models = eval_fitDiscrete_models(phy_new, dat_new, numcol = 1)
#TODO: 2) plot rates of discrete characters (see: http://blog.phytools.org/2016/06/plot-method-for-geigers-fitdiscrete.html)

best_model = names(models[[1]])
dat_vector = models[[2]]

plot_simmmap_pdf(phy_new, dat_vector, best_model, nsim = 100)
