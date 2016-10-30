from genBank import *

species = 'Apteronotidae'
email = 'victor_tagliacollo@yahoo.com.br'
gene  = 'COI'
program = 'muscle'

x = provides_species_get_align(query, gene, email,  batch=500, program = 'muscle')
