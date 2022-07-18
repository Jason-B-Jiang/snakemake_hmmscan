# -----------------------------------------------------------------------------
#
# Generate text file of hmmscan commands for making Pfam domain assignments to
# orthologs in orthogroups.
#
# Jason Jiang - Created: 2022/07/18
#               Last edited: 2022/07/18
#
# Reinke Lab - Microsporidia Orthologs
#
# -----------------------------------------------------------------------------

import sys
import glob
import os
from typing import List

################################################################################

def main() -> None:
    """Creates a text file with all hmmscan commands to run to produce Pfam domain
    assignments for all orthogroups.
    
    Command line arguments:
        $1 = filepath to folder with orthogroup fasta files
        $2 = filepath to file we want to write hmmscan commands to
        $3 = filepath to Pfam HMM library
    """
    og_seqs, hmmscan_cmds, pfam_lib = sys.argv[1:]

    commands = make_hmmscan_command_list(og_seqs, hmmscan_cmds, pfam_lib)

    with open(hmmscan_cmds, 'w') as f:
        f.writelines(commands)

################################################################################

## Helper functions

def make_hmmscan_command_list(og_seqs: str, pfam_lib: str) ->\
    List[str]:
    """Returns a list of hmmscan commands to run, to produce Pfam domain
    assignments for orthogroups.

    Input:
        og_seqs: directory with fasta files for each orthogroup, from Orthofinder
        pfam_lib: filepath to Pfam HMM library (.hmm file)
    """
    commands = []

    for file in glob.glob(f"{og_seqs}/*.fa"):
        og = os.path.basename(file)[:-3]  # get orthogroup name from file
        commands.append(
            f"hmmscan -o ../results/hmmscan/{og} --cut_ga {pfam_lib} {file}\n"
        )
    
    return commands

################################################################################

if __name__ == '__main__':
    main()