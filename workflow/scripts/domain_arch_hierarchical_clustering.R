# -----------------------------------------------------------------------------
#
# Hierarchical clustering of yeast genes by domain architecture conservation
# with orthologs in other species
#
# Jason Jiang - Created: 2022/08/15
#               Last edited: 2022/08/15
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Create a tileplot of domain architecture conservation of yeast genes
#       and their orthologs from other species
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

require(tidyverse)
require(cluster)
require(ggdendro)
require(svglite)

################################################################################

## Global variables

set.seed(42)  # for consistent clustering results

HIERARCHICAL_CLUSTERS_DATA <-
  '../../resources/yeast_genes_hierarchical_clustering.csv'

################################################################################

main <- function() {
  # ----------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = filepath to csv of domain architecture alignments for all single-copy
  #   orthologs to yeast
  #
  #   $2 = filepath to save resulting tileplot image
  # ----------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  aligned_domain_archs <- read_csv(args[1], show_col_types = F)
  out <- args[2]
  
  if (file.exists(HIERARCHICAL_CLUSTERS_DATA)) {
    ortholog_domain_arch_conservation <-
      read_csv(HIERARCHICAL_CLUSTERS_DATA, show_col_types = F)
  } else {
    ortholog_domain_arch_conservation <-
      get_ortholog_domain_arch_conservation(aligned_domain_archs)
    
    write_csv(ortholog_domain_arch_conservation, HIERARCHICAL_CLUSTERS_DATA)
  }

  tileplot <- create_ortholog_tileplot(ortholog_domain_arch_conservation)
  ggsave(plot = tileplot, filename = out, dpi = 600,
         units = 'in', height = 8.49, width = 18.21)
}

################################################################################

## Helper functions

get_ortholog_domain_arch_conservation <- function(aligned_domain_archs) {
  # ----------------------------------------------------------------------------
  # Formats a dataframe of aligned domain architectures for making a tileplot of
  # domain architecture conservation with different yeast orthologs.
  #
  # aligned_DA: dataframe of aligned domain architectures of all microsporidia
  #             species + outgroup species with their yeast orthologs
  # ----------------------------------------------------------------------------
  ortholog_domain_arch_conservation <- aligned_domain_archs %>%
    # group species alphabetically by species names, clades amd yeast orthologs
    mutate(domain_arch_conservation = get_DA_conservation_class(lost_doms,
                                                                gained_doms,
                                                                swapped_doms)) %>%
    select(species, yeast_ortholog, domain_arch_conservation)
  
  # add rows for missing orthologs in each species
  return(fill_missing_orthos(ortholog_domain_arch_conservation))
}


get_DA_conservation_class <- Vectorize(function(lost, gained, swapped) {
  # ----------------------------------------------------------------------------
  # Classify an ortholog to a yeast gene as having identical/conserved
  # domain architectures, lost domains, gained domains or swapped domains.
  # 
  # Domain loss takes precedence over all domain swaps, where an ortholog with
  # domain loss + swap is simply classified as having lost domains only
  # (for plotting purposes)
  #
  # Args:
  #   lost: bool indicating if ortholog has lost domain(s) from yeast
  #
  #   gained: bool indicating if ortholog has gained domain(s) from yeast
  #
  #   swapped: bool indicating if ortholog has swapped domain(s) from yeast
  #
  # ----------------------------------------------------------------------------
  # perfectly conserved DA
  if (all(is.na(c(lost, gained, swapped)))) {
    return('Conserved domain architecture')
  } else if (!is.na(lost)) {
    return('Lost domain(s)')
  } else {
    return('Gained/swapped domain(s)')
  }
}, vectorize.args = c('lost', 'gained', 'swapped'))


fill_missing_orthos <- function(ortholog_domain_arch_conservation) {
  # ----------------------------------------------------------------------------
  # For each species, add rows for yeast genes that the species doesn't have
  # orthologs to, annotating them as "Ortholog not conserved" in the
  # ortholog_domain_arch_conservation dataframe.
  #
  # Args:
  #   ortholog_domain_arch_conservation: dataframe of species + yeast ortholog
  #   domain architecture conservations, with columns for species name, yeast
  #   gene and ortholog domain arch conservation with that yeast gene in the
  #   species.
  #
  # ----------------------------------------------------------------------------
  # get all unique orthologs that microsporidia share with yeast
  yeast_orthos <- unique(ortholog_domain_arch_conservation$yeast_ortholog)
  
  # separate lost_clans_list dataframe into separate dataframes for each species
  sp_DA_conservation <- ortholog_domain_arch_conservation %>%
    group_by(species) %>%
    group_split()
  
  for (i in 1 : length(sp_DA_conservation)) {
    # get orthos that this microsporidia doesn't share with yeast
    missing_orthos <- setdiff(yeast_orthos, sp_DA_conservation[[i]]$yeast_ortholog)
    
    # create dataframe of all missing domains for this species, with a count
    # of zero
    missing_orthos <- data.frame(
      species = unique(sp_DA_conservation[[i]]$species),
      yeast_ortholog = missing_orthos,
      domain_arch_conservation = 'Ortholog not conserved'
    )
    
    # add these missing domains to the original species dataframe
    sp_DA_conservation[[i]] <- rbind(sp_DA_conservation[[i]], missing_orthos)
  }
  
  # join all the split dataframes back together, now with rows added for missing
  # clans in each species
  return(
    bind_rows(sp_DA_conservation) %>%
    arrange(species, yeast_ortholog)
  )
}


create_ortholog_tileplot <- function(ortholog_domain_arch_conservation) {
  sp_order <- c('A_alge', 'A_locu', 'E_aedi', 'E_bien', 'E_brev', 'E_canc',
                'E_cuni', 'E_hell', 'E_hepa', 'E_inte', 'E_roma', 'H_erio',
                'N_ausu', 'N_cera', 'N_cide', 'N_disp', 'N_ferr', 'N_gran',
                'N_homo', 'N_iron', 'N_majo', 'N_mino', 'N_okwa', 'N_pari',
                'O_coll', 'P_epip', 'P_phil', 'S_loph', 'T_cont', 'T_homi',
                'T_rati', 'V_corn', 'V_culi', 'A_sp.', 'M_incu', 'M_daph',
                'P_sacc', 'R_allo', 'H_sapi', 'C_eleg', 'D_reri', 'D_mela',
                'D_disc', 'S_pomb')
  
  agnes_cluster_order <- get_yeast_genes_clusters(ortholog_domain_arch_conservation)
  ortholog_domain_arch_conservation$yeast_ortholog <-
    factor(x = ortholog_domain_arch_conservation$yeast_ortholog,
           levels = ortholog_domain_arch_conservation$yeast_ortholog[agnes_cluster_order],
           ordered = T)
  
  tileplot <- ggplot(ortholog_domain_arch_conservation, aes(x = yeast_ortholog,
                                                            y = species)) +
    geom_tile(aes(fill = domain_arch_conservation)) +
    scale_fill_manual(values=c("cadetblue", "yellow", "red", "black")) +
    scale_y_discrete(limits = rev(sp_order)) +
    labs(x = 'Saccharomyces cerevisiae ortholog',
         title = str_c(length(unique(ortholog_domain_arch_conservation$yeast_ortholog)),
                       ' yeast genes with SCOs in ≥1 other species')) +
    theme(axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 14, color = 'black'),
          axis.title = element_text(size = 18, color = 'black'),
          legend.text = element_text(size = 18, color = 'black'),
          legend.title = element_blank(),
          plot.title = element_text(size = 18, face = 'bold'))
  
  return(tileplot)
}


get_yeast_genes_clusters <- function(ortholog_domain_arch_conservation) {
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # 'widen' the dataframe by making a column for each species, so we have a
  # table of yeast genes and their domain arch conservation with single-copy
  # orthologs from each species
  spread_domain_arch <- ortholog_domain_arch_conservation %>%
    spread(species, domain_arch_conservation)
  
  # convert all character columns to factors, for calculating gower distance
  spread_domain_arch[sapply(spread_domain_arch, is.character)] <-
    lapply(spread_domain_arch[sapply(spread_domain_arch, is.character)], as.factor)
  
  # calculate pairwise similarities between all yeast genes based on their
  # domain architecture conservation status across all other species, using
  # gower's distance (distance metric for mixed data types)
  #
  # shape of distance matrix should be {num yeast genes} x {num yeast genes}
  gower_dist <- daisy(spread_domain_arch, metric = 'gower')
  
  # agglomerative clustering, using ward's linkage for merging clusters
  agnes_clust <- agnes(as.matrix(gower_dist), diss = TRUE, keep.diss = TRUE,
                       method = 'ward')
  
  # dendrogrames
  agnes_dendro <- as.dendrogram(agnes_clust)
  
  # get dendrogram orders
  return(order.dendrogram(agnes_dendro))
}

################################################################################

main()