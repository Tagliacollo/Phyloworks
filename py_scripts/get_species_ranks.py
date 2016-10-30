def get_tax_id(species):
	"""to get data from ncbi taxomomy, we need to have the taxid. we can
    get that by passing the species name to esearch, which will return
    the tax id"""

	# http://stackoverflow.com/questions/16504238/attempting-to-obtain-taxonomic-information-from-biopython
	
	search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
	record = Entrez.read(search)
	search.close()

	return (record['IdList'][0])


def get_tax_data(taxid):
	"""once we have the taxid, we can fetch the record"""

	# http://stackoverflow.com/questions/16504238/attempting-to-obtain-taxonomic-information-from-biopython

	search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
	record = Entrez.read(search)
	search.close()

	return (record)


def get_phylum_name(species):
	# http://stackoverflow.com/questions/16504238/attempting-to-obtain-taxonomic-information-from-biopython
	# (most of it comes from this website)

	sp_rank = defaultdict(list)

	taxid = get_tax_id(species)
	data = get_tax_data(taxid)
	lineage = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank'] in ['order']}

	sp_rank[species].append(lineage)

	return(sp_rank)