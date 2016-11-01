from Bio import Entrez
from Bio import SeqIO, AlignIO
from collections import defaultdict

def provides_species_get_unalign_fasta(species, gene, email, outfilename, batch, rank):

	efetch      = gb_efetch(species, gene, email, batch)
	seqs_dict   = gb_efetch_parse(efetch, rank)
	

	for key, value in seqs_dict.items():
		seqs_dict_by_genus = defaultdict(list)
		seqs_dict_by_genus[key] = value

		fas_unalign = fasta_to_be_aln(seqs_dict_by_genus, gene, outfilename)

	return(fas_unalign)

def gb_efetch(species, gene, email, batch):

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

	print('...done\n')

	return(handle_efetch)

def gb_efetch_parse(efetch_to_be_parsed, rank):
	
	elements = defaultdict(list)

	for seqs in SeqIO.parse(efetch_to_be_parsed, 'gb'):
		key = seqs.annotations['taxonomy'][rank] # change number [-1] if we
		value = (seqs.annotations['organism'], # want seqs for higher rank level (defaut = 'genus') 
				 seqs)

		elements[key].append(value)

	return(elements)

def fasta_to_be_aln(seqs_dict, gene, outfilename):

	seq_list = []
	for key, value in seqs_dict.items():
		for i, seq in enumerate(value):
			seq_list.append(value[i][1])

	SeqIO.write(seq_list, '%s_%s.fasta' % (outfilename, gene), 'fasta')
	
	return (seq_list)
