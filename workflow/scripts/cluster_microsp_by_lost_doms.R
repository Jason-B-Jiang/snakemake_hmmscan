# -----------------------------------------------------------------------------
#
# Cluster microsporidia by domain losses with yeast orthologs
#
# Jason Jiang - Created: 2022/08/15
#               Last edited: 2022/08/15
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Use k-medoids clustering to cluster microsporidia species by their
#       domain losses w.r.t to their yeast orthologs
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(cluster))
suppressMessages(library(Rtsne))

################################################################################

set.seed(42)  # for consistent results

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = filepath to aligned ortholog domain architectures
  #   $2 = filepath to microsporidia + outgroup species clades
  #   $3 = hashmap mapping pfam domains back to their clans
  #   $4 = filepath to save resulting t-SNE plot in
  # ---------------------------------------------------------------------------
  
}

################################################################################

get_species_lost_clans <- function(aligned_domain_archs, fam_to_clan) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    aligned_domain_archs %>%
      select(species, lost_doms) %>%
      filter(!is.na(lost_doms)) %>%
      separate_rows(lost_doms, sep = '; ') %>%
      rowwise() %>%
      mutate(clans = fam_to_clan[[lost_doms]]) %>%
      ungroup() %>%
      select(-lost_doms) %>%
      distinct(.keep_all = T) %>%
      mutate(clan_lost = 1)
  )
}