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


	path ='/Users/Tagliacollo/anaconda/bin/muscle'

	cline = [path, '-quiet', '-maxiters', '1', '-diags']

	#cline = ['muscle', '-quiet', '-maxiters', '1', '-diags']

	child = sub.Popen(cline, stdin  = sub.PIPE, 
							 stdout = sub.PIPE, 
							 stderr = sub.PIPE,
							 universal_newlines = True, 
							 shell = (sys.platform!="win32"))
	
	stdout, stderr = child.communicate()
	
	SeqIO.write(seq_list, child.stdin, "fasta")
	
	child.stdin.close()

	align = AlignIO.read(child.stdout, "clustal")

	return(align)



