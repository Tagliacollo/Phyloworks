from genBank import *
from utilities import output_path
import os

################

email = 'victor_tagliacollo@yahoo.com.br'
gene  = 'COI'
cwd = '/Volumes/VATagliacollo/GitHub/Phyloworks/py_scripts'


ssp_list = ['Gymnotus carapo']

################ FASTA

for sp in ssp_list:
	os.chdir(output_path(cwd, sp))
	print ('acquiring %s sequences for: %s' % (gene, sp))

	provides_species_get_unalign_fasta(sp, gene, email, outfilename=sp, batch=500, rank = -4)


################ ALIGNMENT



################ PHYLOGENY

