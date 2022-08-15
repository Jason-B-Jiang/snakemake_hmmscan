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
  args <- commandArgs(trailingOnly = T)
  aligned_domain_archs <- read_csv(args[1], show_col_types = F)
  clades <- make_clades_hash(args[2])
  fam_to_clan <- read_rds(args[3])
  
  species_clan_losses <- get_species_lost_clans(aligned_domain_archs, fam_to_clan)
  
}

################################################################################

get_species_lost_clans <- function(aligned_domain_archs, fam_to_clan) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    fill_missing_clans(
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
    ) %>%
      pivot_wider(names_from = clans, values_from = clan_lost)
  )
}


fill_missing_clans <- function(clan_losses) {
  # ----------------------------------------------------------------------------
  # For each species, add rows of zero for clans that are not lost or not found
  # in the species' orthologs to yeast
  # ----------------------------------------------------------------------------
  # get all unique orthologs that microsporidia share with yeast
  lost_clans <- unique(clan_losses$clans)
  
  # separate lost_clans_list dataframe into separate dataframes for each species
  sp_clan_losses <- clan_losses %>%
    group_by(species) %>%
    group_split()
  
  for (i in 1 : length(sp_clan_losses)) {
    # get orthos that this microsporidia doesn't share with yeast
    missing_clans <- setdiff(lost_clans, sp_clan_losses[[i]]$clans)
    
    # create dataframe of all missing domains for this species, with a count
    # of zero
    missing_clans <- data.frame(
      species = sp_clan_losses[[i]]$species[1],
      clans = missing_clans,
      clan_lost = 0
    )
    
    # add these missing domains to the original species dataframe
    sp_clan_losses[[i]] <- rbind(sp_clan_losses[[i]], missing_clans)
  }
  
  # join all the split dataframes back together, now with rows added for missing
  # clans in each species
  return(bind_rows(sp_clan_losses))
}