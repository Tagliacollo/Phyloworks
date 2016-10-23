from Bio import Entrez
import pandas as pd

def gb_efetch(ssp_query, email):

	Entrez.email = email

	# Parte 1: egquery
	handle_egquery = Entrez.egquery(term = ssp_query)
	egquery_read   = Entrez.read(handle_egquery)

	for row in egquery_read["eGQueryResult"]:
		if row["DbName"]=="nuccore":
			num_retmax = row["Count"]

	handle_egquery.close()

	# Parte 2: esearch 
	handle_esearch = Entrez.esearch(db = "nucleotide", term = ssp_query, 
									retmax = num_retmax, usehistory = "y")
	esearch_read = Entrez.read(handle_esearch)
	id_list = esearch_read["IdList"]

	handle_esearch.close()

	# Parte 3: epost
	handle_epost = Entrez.epost(db = "nucleotide", id=",".join(map(str,id_list)))
	epost_read = Entrez.read(handle_epost)

	handle_epost.close()

	try:
		handle_efetch = Entrez.efetch(db = "nucleotide", retmode="xml", webenv = epost_read["WebEnv"], 
								 	  query_key = epost_read["QueryKey"])

	except HTTPError as err:
		attempt = 1
		if 500 <= err.code <= 599:
			print("Received error from server %s" % err)
		else:
			raise

	return(handle_efetch)

def gb_efetch_to_df(ssp_query, handle_efetch):

	efetch_read = Entrez.read(handle_efetch)
	
	df = pd.DataFrame()
	for single in efetch_read:
		df = pd.concat([df, gb_efetch_info(ssp_query, single)])

	return(df)

def gb_efetch_info(ssp_query, single_efetch):

	try:
		gene_name   = single_efetch["GBSeq_feature-table"][1]['GBFeature_quals'][0]['GBQualifier_value']
	except:
		gene_name = 'unavailable'
		pass

	try:
		seq = single_efetch['GBSeq_sequence']
	except:
		seq = 'NaN'

	out_name	= ssp_query	
	org_name    = single_efetch["GBSeq_organism"]
	gb_number   = single_efetch["GBSeq_primary-accession"]
	seq_size    = single_efetch['GBSeq_feature-table'][0]['GBFeature_intervals'][0]['GBInterval_to']
	upload_date = single_efetch['GBSeq_update-date']
	paper_title = single_efetch['GBSeq_references'][0]['GBReference_title']
	nucl_seq    = '>%s_%s_%s %s' % (org_name.replace(' ', '_'), gene_name, gb_number, seq)


	dat_list  = [ gene_name,   org_name,   gb_number,   seq_size,  nucl_seq,   upload_date,   paper_title]	
	col_names = ['gene_name', 'org_name', 'gb_number', 'seq_size', 'nucl_seq', 'upload_date', 'paper_title']

	df = pd.DataFrame(data=[dat_list], columns = col_names, index = [out_name])

	return(df)

x = gb_efetch('apteronotus', 'victor_tagliacollo@yahoo.com.br')
y = gb_efetch_to_df('apteronotus', x)


