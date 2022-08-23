# -----------------------------------------------------------------------------
#
# Create hashtable of single-copy ortholog sequence alignments
#
# Jason Jiang - Created: 2022/08/23
#               Last edited: 2022/08/23
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Make a hashtable mapping each protein in an ortholog to an alignment
# with its yeast protein ortholog
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

require(tidyverse)
require(seqinr)
require(glue)
require(Biostrings)

################################################################################

## Global variables

# Parameters for Needleman-Wunsch alignment
GAP_OPEN = 10.0
GAP_EXTEND = 0.5
ALIGNMENT_TYPE = 'global'

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = folder w/ fasta sequences for all orthogroups
  #   $2 = all single-copy ortholog pairs with yeast from orthofinder
  #   $3 = output filepath
  # ---------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  orthogroup_seqs <- args[1]
  orthogroups <- read_csv(args[2], show_col_types = F)
  out <- args[3]
  
  # make hashtable mapping all orthologs in orthogroups to their sequences
  orthogroup_seqs_hash <- make_orthogroup_seqs_hash(orthogroup_seqs,
                                                    unique(orthogroups$orthogroup))
  
  # make hashtable mapping all single-copy orthologs to yeast genes to their
  # alignments with their yeast orthologs
  aligned_orthologs_hash <- make_aligned_orthologs_hash(orthogroups,
                                                        orthogroup_seqs_hash)
  
  # save the alignments hashtable to the specified filepath
  write_rds(aligned_orthologs_hash, out)
}
  
################################################################################

## Helper functions

make_orthogroup_seqs_hash <- function(orthogroup_seqs,
                                      orthogroups_of_interest) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  base_dir <- dirname(list.files(orthogroup_seqs, full.names = T)[1])
  seq_fasta <- unname(sapply(orthogroups_of_interest,
                             function(og) {str_c(base_dir, '/', og, '.fa')}))
  
  # TODO - vectorize this loop
  orthogroup_seqs_hash <- new.env()
  for (seq in seq_fasta) {
    orthogroup_seqs <- new.env()
    fa <- read.fasta(seq, seqtype = 'AA', as.string = T)
    orthologs <- getName(fa)
    
    for (ortholog in orthologs) {
      orthogroup_seqs[[ortholog]] <- fa[[ortholog]][1]
    }
    
    orthogroup_seqs_hash[[basename(str_remove(seq, '\\.fa'))]] <-
      orthogroup_seqs
  }
  
  return(orthogroup_seqs_hash)
}


make_aligned_orthologs_hash <- function(orthogroups, orthogroup_seqs_hash) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  aligned_orthologs_hash <- new.env()
  for (i in 1 : nrow(orthogroups)) {
    aligned_orthologs_hash[[orthogroups$species_ortholog[i]]] <-
      NW_align_orthologs(
        orthogroup_seqs_hash[[orthogroups$orthogroup[i]]][[orthogroups$species_ortholog[i]]],
        orthogroup_seqs_hash[[orthogroups$orthogroup[i]]][[orthogroups$yeast_ortholog[i]]],
        orthogroups$is_microsp[i]
      )
  }
  
  return(aligned_orthologs_hash)
}


NW_align_orthologs <- function(species_seq, yeast_seq, is_microsp) {
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  if (is_microsp) {
    MATRIX = 'BLOSUM45'
  } else {
    MATRIX = 'BLOSUM62'
  }
  
  NW_alignment <- pairwiseAlignment(species_seq, yeast_seq,
                                    type = ALIGNMENT_TYPE,
                                    substitutionMatrix = MATRIX,
                                    gapOpening = GAP_OPEN,
                                    gapExtension = GAP_EXTEND)
  
  return(list(species_alignment = toString(pattern(NW_alignment)),
              yeast_alignment = toString(subject(NW_alignment)),
              percent_id = pid(NW_alignment, type = 'PID1')))
}

################################################################################

main()