query = 'Apteronotidae'
email = 'victor_tagliacollo@yahoo.com.br'

efetch = gb_efetch(query, 'co1', email)
seqs   = gb_efetch_parse(efetch)


for x in seqs.items():
	print(x[0])

