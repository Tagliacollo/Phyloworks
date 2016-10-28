import os
from Bio.Nexus import Nexus

def check_taxa(matrices):
    '''Checks that nexus instances in a list [(name, instance)...] have
        the same taxa, provides useful error if not and returns None if
        everything matches
        From: http://biopython.org/wiki/Concatenate_nexus
    '''
    first_taxa = matrices[0][1].taxlabels
    for name, matrix in matrices[1:]:
        first_only = [t for t in first_taxa if t not in matrix.taxlabels]
        new_only = [t for t in matrix.taxlabels if t not in first_taxa]

        if first_only:
            missing = ', '.join(first_only)
            msg = '%s taxa %s not in martix %s' % (matrices[0][0], missing, name)
            raise Nexus.NexusError(msg)

        elif new_only:
            missing = ', '.join(new_only)
            msg = '%s taxa %s not in all matrices'  % (name, missing)
            raise Nexus.NexusError(msg)

    return None # will only get here if it hasn't thrown an exception


def concat(mypath, same_taxa):
    ''' Combine multiple nexus data matrices in one partitioned file.
        By default this will only work if the same taxa are present in each file
        use  same_taxa=False if you are not concerned by this
        From: http://biopython.org/wiki/Concatenate_nexus
        small change: added onlyfiles block to remove hidden files
    '''

    onlyfiles = []
    for item in os.listdir(mypath):
        if not item.startswith('.') and os.path.isfile(os.path.join(mypath, item)):
            onlyfiles.append(item)
    
    nexi = []
    for nex in onlyfiles:
        nex_open = open(nex, 'r')
        nex_save = Nexus.Nexus(nex_open)
        nexi.append((nex, nex_save))
         
    if same_taxa:
        if not check_taxa(nexi):
            return Nexus.combine(nexi)
    else:
        return Nexus.combine(nexi)


def output_conc_nex(mypath, outfilename, same_taxa=False):

    os.chdir(mypath)
    combined = concat(mypath, same_taxa)
    combined.write_nexus_data(filename=open('%s.nex' % (outfilename), 'w'))

    return None

### how to use it (below)

mypath = '/Users/Tagliacollo/Downloads/sate-gblocks-clean-min-52-taxa'
outfilename = 'Harrington_2016'

output_conc_nex(mypath, outfilename)

