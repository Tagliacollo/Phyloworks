library(phytools)
library(geiger)

# source functions
source('/Volumes/VATagliacollo/GitHub/Phyloworks/r_scripts/R_plot_functions.R')
source('/Volumes/VATagliacollo/GitHub/Phyloworks/r_scripts/R_functions_for_phylo_studies.R')

# set wd
setwd("/Volumes/VATagliacollo/GitHub/Phyloworks")

# load data
phy = read.tree('/Volumes/VATagliacollo/GitHub/Phyloworks/phylo/random_tree.nwk')

dat.csv   = read.csv('/Volumes/VATagliacollo/GitHub/Phyloworks/raw_data/random_Phylomorphospace.csv',
                   row.names = 1)
areas.csv = read.csv('/Volumes/VATagliacollo/GitHub/Phyloworks/raw_data/random_BiogeoArea.csv',
                     row.names = 1)

# clean up data Phylomorphospace 
clean = pw_Cleanup_data(phy, dat.csv)
phy = clean[[1]]
dat.csv = clean[[2]]

# plot phylomorphospace
output = '/Volumes/VATagliacollo/GitHub/Phyloworks/plots/random_phylomorphoscape'
pw_plot_Phylomorphospace(phy, dat.csv, anc.estimates=FALSE, pdf.name=output)


# clean up data BiogeoAreas 
clean = pw_Cleanup_data(phy, areas.csv)
phy = clean[[1]]
dat.csv = clean[[2]]

# plot BiogeoAreas
output = '/Volumes/VATagliacollo/GitHub/Phyloworks/plots/random_BiogeoArea'
pw_plot_BiogeoAreas(phy, areas.csv, output)

