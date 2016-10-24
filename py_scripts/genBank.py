from Bio import Entrez
import pandas as pd
import time

def gb_efetch(taxa, email, batch = 200, gene = 'None'):
	
	Entrez.email = email

	# Parte 1: esearch 
	query = '%s[Orgn]' % (taxa)

	if gene != 'None':
		query = '%s[Orgn] AND %s[Gene]' % (taxa, gene)

	handle_esearch = Entrez.esearch(db = 'nucleotide', term = query, retmax = 100000000)
	esearch_read = Entrez.read(handle_esearch)
	
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
			handle_efetch = Entrez.efetch(db = "nucleotide", retmode = 'xml', rettype = 'gb', 
										  webenv = epost_read["WebEnv"], query_key = epost_read["QueryKey"])

		except HTTPError as err:
			if 500 <= err.code <= 599:
				print("Received error from server %s" % err)
			else:
				raise

	return(handle_efetch)

def gb_efetch_to_df(taxa, handle_efetch):

	efetch_read = Entrez.read(handle_efetch)
	
	df = pd.DataFrame()
	for single in efetch_read:
		df = pd.concat([df, gb_efetch_parse(taxa, single)])

	return(df)

def gb_efetch_parse(ssp_query, single_efetch):

	try:
		gene_name   = single_efetch["GBSeq_feature-table"][1]['GBFeature_quals'][0]['GBQualifier_value']
	except:
		gene_name = 'NaN'
		pass

	try:
		seq = single_efetch['GBSeq_sequence']
	except:
		seq = 'NaN'
		pass

	org_name    = single_efetch["GBSeq_organism"]
	gb_number   = single_efetch["GBSeq_primary-accession"]
	seq_size    = single_efetch['GBSeq_feature-table'][0]['GBFeature_intervals'][0]['GBInterval_to']
	upload_date = single_efetch['GBSeq_update-date']
	paper_title = single_efetch['GBSeq_references'][0]['GBReference_title']
	nucl_seq    = '>%s_%s_%s %s' % (org_name.replace(' ', '_'), gene_name, gb_number, seq)


	dat_list  = [ gene_name,   org_name,   gb_number,   seq_size,  nucl_seq,   upload_date,   paper_title]	
	col_names = ['gene_name', 'org_name', 'gb_number', 'seq_size', 'nucl_seq', 'upload_date', 'paper_title']

	df = pd.DataFrame(data=[dat_list], columns = col_names, index = [gene_name])

	return(df)




ssp_query = 'Gymnotus'
email     = 'victor_tagliacollo@yahoo.com.br'

y = gb_efetch(ssp_query, email, gene = 'COI')
z = gb_efetch_to_df(ssp_query, y)


