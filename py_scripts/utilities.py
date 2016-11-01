import os
from pathlib2 import Path


def output_path(cwd, folder_name):

    repo_dir = Path(cwd).parents[0]
    seqs_dir = os.path.join(str(repo_dir), 
                        "projects/VATagliacollo/phylo_of_gymn/sequences")
    output_path = os.path.join(seqs_dir, folder_name)
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    return (output_path)