import Bio.Align
from Bio.Nexus import Nexus
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
from io import StringIO


def processed_data(aln_nex, cut, outfilename):

    alns = charset_aln(aln_nex)
    
    ls_alns = []
    for aln in alns:
        gaps = cut_gap_in_blocks(aln[1], cut)
        if gaps.get_alignment_length() == 0:
            continue
        else:
            dupl = remove_dup_seqs(gaps)
            ls_alns.append((aln[0], dupl))

    export_nexus(ls_alns, outfilename)

def charset_aln(aln_nex):
# aln: 	nexus aln with charasets
# return: a list of generic biopython alns, each representing an charset
	
	dat = Nexus.Nexus()
	dat.read(aln_nex)
	aln_nex = AlignIO.read(open(aln_nex), "nexus")
    
	aln_by_charset   = []
	for name in dat.charsets:
		sites = dat.charsets[name]
		start = min(sites)
		stop  = max(sites) + 1
		aln_by_charset.append((name, aln_nex[:, start:stop]))

	return(aln_by_charset)

def cut_gap_in_blocks(aln, cut):
# From - https://stackoverflow.com/questions/28220301/python-remove-special-column-from-multiple-sequence-alignment
    n = float(len(aln[0]))
    i = 0
    while i < n:
        ct = 0
        while i+ct < n and aln[:, i+ct].count('?') / float(len(aln[:, i+ct])) > cut:
            ct += 1

        if ct > 0:     #  delete columns [i:i+ct]
            if i == 0:
                aln = aln[:, ct:]
            elif i+ct == n:
                aln = aln[:, :i]
            else:
                aln = aln[:, :i] + aln[:, i+ct:]
            n -= ct    #  seq. ct positions shorter
        else:          # nothing to delete, proceed
            i += 1

    return(aln)


def remove_dup_seqs(aln):
# From - http://lists.open-bio.org/pipermail/biopython/2010-April/012615.html  

    aln_no_dup = []
    checksums  = []
    for record in aln:
        checksum = seguid(record.seq)
        if checksum in checksums:
            #print ("Ignoring %s" % record.id)
            continue
        else:
            checksums.append(checksum)
            aln_no_dup.append(record)

    return(Bio.Align.MultipleSeqAlignment(aln_no_dup))


def export_nexus(ls_alns, outfilename):

    nexi = []
    for aln in ls_alns:
        output = StringIO()
        AlignIO.write(aln[1], output, 'nexus')
        p = Nexus.Nexus()
        p.read(output.getvalue())
    
        nexi.append((aln[0], p))

    combined = Nexus.combine(nexi)
    combined.write_nexus_data(outfilename, interleave=True)



########## WORKING

data = '/Users/Tagliacollo/Google-Drive/05-G-forms/raw-data/aln/mtx/GymnSuperMTXAll.nex'
#data = '/Users/Tagliacollo/Google-Drive/05-G-forms/raw-data/aln/mtx/TestMTX.nex'
cut  = 0.60
outfilename = '/Users/Tagliacollo/Google-Drive/05-G-forms/raw-data/aln/mtx/processed/GymnSuperMTXAll60.nex'

processed_data(data, cut, outfilename)


