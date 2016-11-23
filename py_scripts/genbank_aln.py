from os import listdir
from os.path import isfile, join
from Bio.Align.Applications import ClustalOmegaCommandline
import os



path = "/Volumes/VATagliacollo/GitHub/Phyloworks/projects/VATagliacollo/phylo_of_gymn/sequences/Gymnotiformes"

fastas = [f for f in listdir(mypath) if isfile(join(mypath, f))]

os.chdir('/Volumes/VATagliacollo/GitHub/Phyloworks/projects/VATagliacollo/phylo_of_gymn/sequences/Gymnotiformes')

for fasta in fastas:
	print(fasta)
	muscle_cline = MuscleCommandline(input = fas, out= '%s.aln' % (fasta))


def get_clustalo_aln(in_fas, out_fas)







from subprocess import Popen, PIPE, STDOUT
from Bio import SeqIO, AlignIO
import sys

# check this link to install to install argtable 2 and clustal: https://www.biostars.org/p/128261/
# best anwser is given by Alex Reynolds
# **Important** - muscle and mafft must be installed in anaconda python. 
# To install muscle type on the command line: conda install -c tomkinsc muscle=3.8.31
# To install mafft type on the command line: conda install -c bioconda mafft=7.221
# To install clustalOmega on the command line: conda install -c etetoolkit argtable2=2.13

records = (r for r in SeqIO.write("CYTB_Fritz.fas", "fasta") if len(r) < 900)

cline = ['clustalo', '-', '--auto', '-v']

proc = Popen(cline, stdout=PIPE, stdin=PIPE, stderr=STDOUT)

stdout = proc.communicate(input=records)



def align_seq_cmd(seq_list, program = 'muscle'):

	# **Important** - muscle and mafft must be installed in anaconda python. 
	# To install muscle type on the command line: conda install -c tomkinsc muscle=3.8.31
	# To install mafft type on the command line: conda install -c bioconda mafft=7.221
	# To install clustalOmega on the command line: conda install -c etetoolkit argtable2=2.13

	# Muscle arguments
	#   -in <inputfile>    Input file in FASTA format (default stdin)
	#    -out <outputfile>  Output alignment in FASTA format (default stdout)
	#    -diags             Find diagonals (faster for similar sequences)
	#    -maxiters <n>      Maximum number of iterations (integer, default 16)
	#    -maxhours <h>      Maximum time to iterate in hours (default no limit)
	#    -html              Write output in HTML format (default FASTA)
	#    -msf               Write output in GCG MSF format (default FASTA)
	#    -clw               Write output in CLUSTALW format (default FASTA)
	#    -clwstrict         As -clw, with 'CLUSTAL W (1.81)' header
	#    -log[a] <logfile>  Log to file (append if -loga, overwrite if -log)
	#    -quiet             Do not write progress messages to stderr
	#    -version           Display version information and exit

	# Without refinement (very fast, avg accuracy similar to T-Coffee): -maxiters 2
	# Fastest possible (amino acids): -maxiters 1 -diags -sv -distance1 kbit20_3
	# Fastest possible (nucleotides): -maxiters 1 -diags

	handle = StringIO()
	SeqIO.write(seq_list, handle, 'fasta')

	handle.getvalue()

	for record in SeqIO.parse(handle.getvalue(), "fasta"):
		print(record.id)


	x = open('Apteronotus.fasta', 'r')
	x.read()

records = (r for r in SeqIO.parse("CYTB_Fritz.fas", "fasta") if len(r) < 900)

cline   = ['muscle', '-quiet', '-maxiters', '1', '-diags']

#cline = ClustalOmegaCommandline(infile=records, outfile=out_file, verbose=True, auto=True)
#cline = ClustalOmegaCommandline(infile=records, outfile=out_file, verbose=True, auto=True)

#cline = ['clustalo', '--auto', '-v']

child = Popen(cline, stdin = PIPE,
					 stdout = PIPE,
					 stderr = PIPE,
					 universal_newlines = True) 

SeqIO.write(records, child.stdin, "fasta")
align = AlignIO.read(child.stdout, 'fasta')
	
stdin, stdout = child.communicate() 


	child.stdin.close()
	child.stdout.close()


import subprocess
from Bio import AlignIO

from Bio.Align.Applications import MuscleCommandline
muscle_cline = MuscleCommandline(input="CYTB_Fritz.fas")
child = subprocess.Popen(str(muscle_cline),
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          universal_newlines=True,
                          shell=(sys.platform!="win32"))

align = AlignIO.read(child.stdout, "fasta")
print(align)
