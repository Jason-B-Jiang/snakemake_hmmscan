# -----------------------------------------------------------------------------
#
# Clean microsporidia name + host data for predictions
#
# Jason Jiang - Created: 2022/07/29
#               Last edited: 2022/09/09
#
# Reinke Lab - Microsporidia orthologs
#
# Clean up microsporidia/host name data for predicting with rule-based methods,
# and add info about where each microsporidia/host name are found in text.
#
# -----------------------------------------------------------------------------

library(tidyverse)
library(readxl)

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  
  # TODO - replace w/ command line arguments
  sgd_to_uniprot_names <- read_rds('../../data/yeast_genome_to_uniprot.rds')
  ptc_data <- format_ptc_data(readxl::read_xls('../../data/supp_11.xls'),
                              sgd_to_uniprot_names)
  ortholog_alignments <- read_rds('../../resources/ortholog_alignments.rds')
  orthogroups <- read_csv('../../data/orthogroups.csv')
}

################################################################################

## Helper functions

format_ptc_data <- function(ptc_data, sgd_to_uniprot_names) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    ptc_data %>%
      select(Gene, CDS_length, dist_from_CDS_end, `HMM PTC classification`) %>%
      rename(gene = Gene, cds_length = CDS_length,
             dist_from_cds_end = dist_from_CDS_end,
             ptc_tolerance = `HMM PTC classification`) %>%
      filter(!is.na(ptc_tolerance)) %>%
      mutate(ptc_position = cds_length - dist_from_cds_end) %>%
      group_by(gene) %>%
      arrange(ptc_position, .by_group = T) %>%
      mutate(ptc_tolerance = str_c(ptc_tolerance, collapse = ', '),
             ptc_position = str_c(as.character(ptc_position), collapse = ', ')) %>%
      select(-dist_from_cds_end) %>%
      distinct(.keep_all = T) %>%
      ungroup() %>%
      mutate(critical_ptc_position = get_critical_ptc_position(ptc_tolerance,
                                                               ptc_position)) %>%
      rowwise() %>%
      mutate(gene = sgd_to_uniprot_names[[gene]])
  )
}


get_critical_ptc_position <- Vectorize(function(ptc_tolerance, ptc_position) {
  # ---------------------------------------------------------------------------
  # Return corresponding position of the last lethal PTC in a gene
  # ---------------------------------------------------------------------------
  ptc_tolerance <- str_split(ptc_tolerance, ', ')[[1]]
  ptc_position <- str_split(ptc_position, ', ')[[1]]
  
  # corresponding ptc position in gene for last lethal ptc
  critical_ptc_position <- ptc_position[tail(which(ptc_tolerance == 'D'), n = 1)]
  
  if (length(critical_ptc_position) == 0) {
    # all PTCs in this gene are tolerated
    return(NA)
  }
  
  return(critical_ptc_position)
}, vectorize.args = c('ptc_tolerance', 'ptc_position'))


get_lost_residues_in_species_ortholog <- function(species_alignment,
                                                  yeast_alignment) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  i_yeast <- 0  # corresponding yeast residue in the ortholog alignment
  lost_residues <- integer()
  
  for (i in 1 : str_length(species_alignment)) {
    species_aln <- species_alignment[i]
    yeast_aln <- yeast_alignment[i]
    
    if (yeast_aln != '-' & microsp_aln != '-') {
      # yeast residue is aligned to a non-gap position in the species ortholog,
      # so add 1 to update the current yeast residue in the alignment
      i_yeast <- i_yeast + 1
    }
    
    if (microsp_aln == '-' & yeast_aln != '-') {
      lost_residues <- append(lost_residues, i_yeast)
    }
  }
  
  return(ifelse(length(lost_residues) == 0,
                NA,
                str_c(lost_residues, collapse = '; ')))
}

################################################################################