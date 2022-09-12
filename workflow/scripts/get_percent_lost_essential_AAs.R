# -----------------------------------------------------------------------------
#
# Clean microsporidia name + host data for predictions
#
# Jason Jiang - Created: 2022/07/29
#               Last edited: 2022/09/12
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
  
  # add columns to orthogroups dataframe indicating number of lost residues
  # in species orthologs to yeast, and number of such residues that are
  # essential
  orthogroups <- annotate_essential_lost_residues(orthogroups, ptc_data,
                                                  ortholog_alignments)
}

################################################################################

## Helper functions

format_ptc_data <- function(ptc_data, sgd_to_uniprot_names) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    ptc_data %>%
      select(Gene, CDS_length, dist_from_CDS_end, `HMM PTC classification`) %>%
      rename(gene = Gene, protein_length = CDS_length,
             dist_from_cds_end = dist_from_CDS_end,
             ptc_tolerance = `HMM PTC classification`) %>%
      filter(!is.na(ptc_tolerance)) %>%
      mutate(ptc_position = protein_length - dist_from_cds_end,
             # original CDS_length column included stop codon in length
             protein_length = protein_length - 1) %>%
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


annotate_essential_lost_residues <- function(orthogroups, ptc_data,
                                             ortholog_alignments) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    orthogroups %>%
      filter(yeast_ortholog %in% ptc_data$gene) %>%
      select(orthogroup, species, is_microsp, yeast_ortholog, species_ortholog,
             yeast_len, species_len) %>%
      rename(gene = yeast_ortholog) %>%
      # add columns for ptc data of each gene
      left_join(ptc_data, by = 'gene') %>%
      rename(yeast_ortholog = gene) %>%
      rowwise() %>%
      mutate(species_lost_residues =
               get_lost_residues_in_species_ortholog(ortholog_alignments[[species_ortholog]]$species_alignment,
                                                     ortholog_alignments[[species_ortholog]]$yeast_alignment)) %>%
      ungroup() %>%
      mutate()
  )
}


get_lost_residues_in_species_ortholog <- function(species_alignment,
                                                  yeast_alignment) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  species_alignment <- str_split(species_alignment, '')[[1]]
  yeast_alignment <- str_split(yeast_alignment, '')[[1]]
  i_yeast <- 0  # corresponding yeast residue in the ortholog alignment
  lost_residues <- integer()
  
  for (i in 1 : length(species_alignment)) {
    species_aln <- species_alignment[i]
    yeast_aln <- yeast_alignment[i]
    
    if (yeast_aln != '-') {
      # yeast residue is aligned to a non-gap position in the species ortholog,
      # so add 1 to update the current yeast residue in the alignment
      i_yeast <- i_yeast + 1
    }
    
    if (species_aln == '-' & yeast_aln != '-') {
      lost_residues <- append(lost_residues, i_yeast)
    }
  }
  
  return(ifelse(length(lost_residues) == 0,
                NA,
                str_c(lost_residues, collapse = '; ')))
}


classify_lost_amino_acid <- Vectorize(function(lost_residue, ptc_tolerance,
                                     ptc_position, critical_ptc_position) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  if (!str_detect(ptc_tolerance, 'A')) {
    # all lethal PTCs
    if (lost_residue <= critical_ptc_position) {
      return('essential')
    }
    return('unclassified')
  }
  
  if (is.na(critical_ptc_position)) {
    # all PTCs are tolerated
    first_ptc <- as.integer(str_split(ptc_position, '; ')[[1]][1])
    if (lost_residue >= first_ptc) {
      return('dispensable')
    }
    return('unclassified')
  }
  
  # at least one lethal PTC occurs in the gene, and not all PTCs are lethal
  if (lost_residue <= critical_ptc) {
    return('essential')
  }
  
  return('dispensable')
}, vectorize.args = c('lost_residue', 'ptc_tolerance', 'ptc_position',
                      'critical_ptc_position'))

################################################################################