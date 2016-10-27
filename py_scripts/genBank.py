from Bio import Entrez
from Bio import SeqIO, AlignIO
from collections import defaultdict
from Bio.Align.Applications import MuscleCommandline
import subprocess
import sys


def gb_efetch(taxa, gene, email, batch=500):

	try:
		from urllib.error import HTTPError  # for Python 3
	except ImportError:
		from urllib2 import HTTPError  # for Python 2

	Entrez.email = email
	query = '%s[Orgn] AND %s[Gene]' % (taxa, gene)

	#part 1: esearch
	handle_esearch = Entrez.esearch(db = 'nucleotide', term = query, retmax = 100000000)
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
			     seqs.format('fasta'))

		elements[key].append(value)

	return(elements)








