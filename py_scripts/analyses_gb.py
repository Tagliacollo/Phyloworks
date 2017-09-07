from genBank import *
from utilities import output_path
from collections import defaultdict
from Bio import SeqIO
import itertools
import os

################

email = 'victor_tagliacollo@yahoo.com.br'
cwd = '/Volumes/VATagliacollo/GitHub/Phyloworks/py_scripts'


################

ssp_list  = ['Gymnotiformes']
gene_list = [['ND1','MTND1','NADH1'],	# MT
			['ND2','MTND2','NADH2'],	# MT
			['ND3','MTND3','NADH3'],	# MT
			['ND4','MTND4','NADH4'],	# MT
			['ND4L','MTND4L','NADH4L'],	# MT
			['ND5','MTND5','NADH5'],	# MT
			['ND6','MTND6','NADH6'],	# MT
			['CYTB','COB','MTCYB'],		# MT
			['COX1','COI','CO1','COX1','MTCO1','cytochrome oxidase 1'],	# MT
			['CO2','COII','COXII','MTCO2','cytochrome oxidase 2'],		# MT
			['CO3','COIII','COXIII','MTCO3','cytochrome oxidase 3'],	# MT
			['ATP6','ATPASE6','MTATP6','OLI2','OLI4','PHO1'],			# MT
			['ATP8','ATPASE8', 'MTATP8','AAP1'],						# MT
			['ZIC1','ZIC','ZNF201'],				#NUC
			['MYH6','amhc', 'Myhca', 'Myhc-a'],		#NUC
			['RYR3'],								#NUC
			['PTR','ptrA'],							#NUC
			['TBR1'],								#NUC
			['ENC1','NRPB', 'CCL28', 'ENC-1', 'PIG10', 'KLHL35', 'KLHL37', 'TP53I10'],	#NUC
			['GYLT','pomgnt2'],							#NUC
			['SH3PX3','SNX33', 'SNX30', 'SH3PXD3C'],	#NUC
			['PLAGL2','ZNF900'],						#NUC
			['SREB2','GPR85', 'SREB'],					#NUC 
			['RAG1', 'Rag-1'],							#NUC
			['RAG2','Rag-2'],							#NUC
			['RHOD', 'Rho', 'ARHD', 'RHOM', 'RHOHP1'],	#NUC
			['EGR3']]	#NUC


################ FASTA
dict_of_seqRecords = defaultdict(list)

for sp in ssp_list:
	path = ('projects/VATagliacollo/phylo_of_gymn/sequences/%s' % sp)
	os.chdir(output_path(cwd, path))
	for gene in gene_list:
		for symn in gene:
			print ('acquiring %s sequences for: %s' % (symn, sp))
			try:
				dict_of_seqRecords[gene[0]].append(get_seqRecords(sp, symn, email, batch = 500, rank = 1))
			except:
				print('... none found\n')
				pass

	for key, value in dict_of_seqRecords.items():
		if len(value) == 0:
			print('none sequences avalaible for %s ... skipping' % (key))
		else:
			merge_seqs = list(itertools.chain(*value))
			SeqIO.write(merge_seqs, '%s.fas' % key, 'fasta')
	
################ ALIGNMENT


################ PHYLOGENY

