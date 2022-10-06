# -----------------------------------------------------------------------------
#
# Clean microsporidia name + host data for predictions
#
# Jason Jiang - Created: 2022/07/29
#               Last edited: 2022/09/13
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

# Ordering of microsporidia/outgroup species
SP_ORDER <- c('A_alge', 'A_locu', 'E_aedi', 'E_bien', 'E_brev', 'E_canc',
              'E_cuni', 'E_hell', 'E_hepa', 'E_inte', 'E_roma', 'H_erio',
              'N_ausu', 'N_cera', 'N_cide', 'N_disp', 'N_ferr', 'N_gran',
              'N_homo', 'N_iron', 'N_majo', 'N_mino', 'N_okwa', 'N_pari',
              'O_coll', 'P_epip', 'P_phil', 'S_loph', 'T_cont', 'T_homi',
              'T_rati', 'V_corn', 'V_culi', 'A_sp.', 'M_incu', 'M_daph',
              'P_sacc', 'R_allo', 'C_eleg', 'D_disc', 'D_mela', 'D_reri',
              'H_sapi', 'S_pomb')

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
  out <- '../../results/lost_residues_essential.svg'
  
  # add columns to orthogroups dataframe indicating number of lost residues
  # in species orthologs to yeast, and number of such residues that are
  # essential
  orthogroups <- annotate_essential_lost_residues(orthogroups, ptc_data,
                                                  ortholog_alignments)
  
  plot_percent_lost_residues(orthogroups, out)
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
                                                     ortholog_alignments[[species_ortholog]]$yeast_alignment),
             essential_dispensable_unclassified = get_lost_amino_acid_types(species_lost_residues,
                                                                            critical_ptc_position,
                                                                            ptc_tolerance,
                                                                            ptc_position)) %>%
      ungroup() %>%
      separate(essential_dispensable_unclassified,
               c('num_lost_essential', 'num_lost_dispensable', 'num_lost_unclassified'),
               sep = ', ') %>%
      mutate(num_lost = ifelse(!is.na(species_lost_residues),
                                        lengths(str_split(species_lost_residues, ', ')),
                                        0),
             num_lost_essential = as.integer(num_lost_essential),
             num_lost_dispensable = as.integer(num_lost_dispensable),
             num_lost_unclassified = as.integer(num_lost_unclassified)) %>%
      filter(species_len < yeast_len, !is.na(species_lost_residues)) %>%
      # now, summarise average number of amino acid type losses per ortholog
      # for each species
      group_by(species, is_microsp) %>%
      summarise(essential = mean(num_lost_essential),
                dispensable = mean(num_lost_dispensable),
                unclassified = mean(num_lost_unclassified)) %>%
      pivot_longer(cols = c('essential', 'dispensable', 'unclassified'),
                   names_to = 'lost_residue_type',
                   values_to = 'avg_loss_per_ortholog')
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
                str_c(lost_residues, collapse = ', ')))
}


get_lost_amino_acid_types <- function(species_lost_residues,
                                      critical_ptc_position,
                                      ptc_tolerance,
                                      ptc_position) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  if (is.na(species_lost_residues)) {
    return(NA)
  }
  
  species_lost_residues <- as.integer(str_split(species_lost_residues, ', ')[[1]])
  lost_residue_types <- sapply(species_lost_residues,
                               function(a) {classify_lost_amino_acid(a,
                                                                     ptc_tolerance,
                                                                     ptc_position,
                                                                     critical_ptc_position)})
  
  num_essential <- length(lost_residue_types[lost_residue_types == 'essential'])
  num_dispensable <- length(lost_residue_types[lost_residue_types == 'dispensable'])
  num_unclassified <- length(lost_residue_types[lost_residue_types == 'unclassified'])
  
  return(str_c(num_essential, num_dispensable, num_unclassified, sep = ', '))
}


classify_lost_amino_acid <- function(lost_residue, ptc_tolerance,
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
    first_ptc <- as.integer(str_split(ptc_position, ', ')[[1]][1])
    if (lost_residue >= first_ptc) {
      return('dispensable')
    }
    return('unclassified')
  }
  
  # at least one lethal PTC occurs in the gene, and not all PTCs are lethal
  if (lost_residue <= critical_ptc_position) {
    return('essential')
  }
  
  return('dispensable')
}


plot_percent_lost_residues <- function(orthogroups, out) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # Conclusion: Microspordia proteins get shorter (show other fig as supp), but
  # still lose same proportion of essential amino acids
  #
  # Microsporidian proteins likely highly diverged in function to yeast orthologs,
  # relative to outgroup species
  #
  # Most lost residues still essential due to bias for classifying residues as
  # essential
  ggplot(orthogroups, aes(x = species, y = avg_loss_per_ortholog, fill = lost_residue_type)) +
    geom_bar(stat = 'identity') +
    labs(y = 'Average residues lost per ortholog') +
    scale_x_discrete(limits = rev(SP_ORDER)) +
    coord_flip() +
    theme_bw() +
    theme(legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 14, color = 'black'),
          legend.text = element_text(size = 14)) 
}

################################################################################