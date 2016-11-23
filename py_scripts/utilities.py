import os
from pathlib2 import Path
from Bio.SeqUtils.CheckSum import seguid
from Bio import SeqIO


def output_path(cwd, folderpath):

    repo_dir = Path(cwd).parents[0]
    output_path = os.path.join(str(repo_dir), folderpath) 
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    return (output_path)