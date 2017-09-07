# Imports
from Bio import Entrez
from Bio import SeqIO
import io

Entrez.email     = "victor_tagliacollo@yahoo.com.br"  # Always tell NCBI who you are
genomeAccessions = ['NC_015754.1','AP011978.1', 'KX058570.1', 'NC_015078.1', 'AP011979.1', 'NC_004701.1']
search           = " ".join(genomeAccessions)
handle           = Entrez.read(Entrez.esearch(db="nucleotide", term=search, retmode="xml"))
genomeIds        = handle['IdList']
records          = Entrez.efetch(db="nucleotide", id=genomeIds, rettype="gb", retmode="text")

###############################
# Generate Genome Fasta files #
###############################

sequences   = []  # store your sequences in a list
headers     = []  # store genome names in a list (db_xref ids)

for i, record in enumerate(records):

    file_out = open("genBankRecord_"+str(i)+".gb", "w")    # store each genomes .gb in separate files
    file_out.write(record.read())
    file_out.close()

    genomeGenbank   = SeqIO.read("genBankRecord"+str(i)+".gb", "genbank")  # parse in the genbank files
    header         = genome.features[0].qualifiers['db_xref'][0]          # name the genome using db_xfred ID
    sequence       = genome.seq.tostring()                                # obtain genome sequence

    headers.append('>'+header)  # store genome name in list                                     
    sequences.append(sequence)  # store sequence in list

    fasta_out = open("genome"+str(i)+".fasta","w")     # store each genomes .fasta in separate files
    fasta_out.write(header)    # >header ... followed by:
    fasta_out.write(sequence)  # sequence ... 
    fasta_out.close()          # close that .fasta file and move on to next genome
records.close()


# from: http://www.christinabergey.com/notes_code/genbank_to_gene_fastas/gb_to_gene_fastas_CMB.py.txt

from Bio import SeqIO
import os
import sys





for record in SeqIO.parse(handle, 'genbank') :	
	seq = str(record.seq)
	for feature in record.features:
		print feature.type
		print feature.location	
		if feature.type == 'CDS' or feature.type == 'rRNA':
				
			if 'gene' in feature.qualifiers:
				geneName = feature.qualifiers['gene'][0]
			elif 'product' in feature.qualifiers:
				geneName = feature.qualifiers['product'][0]
			else:
				print 'ERROR when parsing feature:'
				print feature.qualifiers
				quit()
			
			geneName = geneName.replace(' ', '-');
			print geneName
			
			geneFile = open(filebase + '_' + record.id + '_'+ geneName + '.fa', 'w')
			geneFile.write('>')
			geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
			geneFile.write(seq[feature.location.start.position:feature.location.end.position])
			geneFile.write("\n\n")
		
		print '----------------------------'





from Bio import Entrez
from Bio import SeqIO
import io

Entrez.email     = "victor_tagliacollo@yahoo.com.br"  # Always tell NCBI who you are
genomeAccessions = ['AP011978.1', 'AP011979.1', 'KX058570.1', 'NC_004692.1', 'NC_004701.1', 'NC_015078.1', 'NC_015754.1']
search           = " ".join(genomeAccessions)
handle           = Entrez.read(Entrez.esearch(db="nucleotide", term=search, retmode="xml"))
genomeIds        = handle['IdList']
records          = Entrez.efetch(db="nucleotide", id=genomeIds, rettype="gb", retmode="text")


for record in SeqIO.parse(records, 'genbank'):
  seq = str(record.seq)
  ID = str(record.id)
  sp = record.description

  for feature in record.features:
  	if feature.type == 'CDS' or feature.type == 'rRNA' or feature.type == 'tRNA' or feature.type == 'D-loop':
  		if 'gene' in feature.qualifiers:
  			geneName = feature.qualifiers['gene'][0]
  		elif 'product' in feature.qualifiers:
  			geneName = feature.qualifiers['product'][0]
  		else:
  			print ('ERROR when parsing feature:')
  			print (feature.qualifiers)
  			quit()

  		geneFile = open(geneName + '.fas', 'a') 
  		geneFile.write('>')	
  		geneFile.write('genome_' + ID + '_' + sp + '_' + geneName + "\n")
  		
  		geneFile.write(seq[feature.location.start.position:feature.location.end.position])
  		geneFile.write("\n\n")



