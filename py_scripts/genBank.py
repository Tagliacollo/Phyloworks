# from get_species_ranks import *
from Bio import Entrez
from Bio import SeqIO, AlignIO
from collections import defaultdict
from Bio.Align.Applications import MuscleCommandline, MafftCommandline, ClustalOmegaCommandline
import subprocess
import sys


def provides_species_get_align(species, gene, email, batch=500, program = 'muscle'):

	efetch      = gb_efetch(species, gene, email, batch=500)
	seqs_dict   = gb_efetch_parse(efetch)
	
	align_multi_seqs  = defaultdict(list)

	for key, value in seqs_dict.items():
		seqs_dict_by_genus = defaultdict(list)
		seqs_dict_by_genus[key] = value

		fas_unalign = fasta_to_be_aln(seqs_dict_by_genus)

		if len(fas_unalign) < 2:
			align_multi_seqs[key].append(fas_unalign)
			
		else:
			align_multi_seqs[key].append(align_seq_cmd(fas_unalign, program))

	return(align_multi_seqs)

def gb_efetch(species, gene, email, batch=500):

	try:
		from urllib.error import HTTPError  # for Python 3
	except ImportError:
		from urllib2 import HTTPError  # for Python 2

	Entrez.email = email
	name = '%s[Orgn] AND %s[Gene]' % (species, gene)

	#part 1: esearch
	handle_esearch = Entrez.esearch(db = 'nucleotide', term = name, retmax = 100000000)
	esearch_read = Entrez.read(handle_esearch)
	handle_esearch.close()
	
	id_list = esearch_read['IdList']
	count   = esearch_read['Count']

	# Parte 2: epost
	handle_epost = Entrez.epost(db = 'nucleotide', id=",".join(map(str,id_list)))
	epost_read = Entrez.read(handle_epost)
	handle_epost.close()

	# Parte 3: efetch
	for start in range(0, int(count), batch):
		end = min(int(count), start + batch)
		print("Going to download record %i to %i" % (start + 1, end))

		try:
			handle_efetch = Entrez.efetch(db = "nucleotide", retmode = 'text', rettype = 'gb', 
										  webenv = epost_read["WebEnv"], query_key = epost_read["QueryKey"])

		except HTTPError as err:
			if 500 <= err.code <= 599:
				print("Received error from server %s" % err)
			else:
				raise

	return(handle_efetch)

def gb_efetch_parse(efetch_to_be_parsed):
	
	elements = defaultdict(list)

	for seqs in SeqIO.parse(efetch_to_be_parsed, 'gb'):
		key = seqs.annotations['taxonomy'][-1]
		value = (seqs.annotations['organism'], 
				 seqs,
				 seqs.annotations['taxonomy'])

		elements[key].append(value)

	return(elements)

def fasta_to_be_aln(seqs_dict):

	seq_list = []
	for key, value in seqs_dict.items():
		for i, seq in enumerate(value):
			seq_list.append(value[i][1])

	return(tuple(seq_list))

def align_seq_cmd(tuple_list, program = 'muscle'):

	# **Important** - muscle and mafft must be installed in anaconda python. 
	# To install muscle type on the command line: conda install -c tomkinsc muscle=3.8.31
	# To install mafft type on the command line: conda install -c bioconda mafft=7.221
	# To install clustalOmega on the command line: conda install -c etetoolkit argtable2=2.13



	if program == 'muscle':
		cline = MuscleCommandline(clwstrict=True)

	if program == 'clustalOmega':
		cline = ClustalOmegaCommandline(verbose=True, auto=True)

	if program == 'mafft':
		cline = MafftCommandline()

	child = subprocess.Popen(str(cline), stdin=subprocess.PIPE, 
							 stdout=subprocess.PIPE, stderr=subprocess.PIPE,
							 universal_newlines=True, shell=(sys.platform!="win32"))

	SeqIO.write(tuple_list, child.stdin, "fasta")
	child.stdin.close()

	align = AlignIO.read(child.stdout, "clustal")

	return(align)