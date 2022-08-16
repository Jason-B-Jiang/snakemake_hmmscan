# -----------------------------------------------------------------------------
#
# Cluster microsporidia by domain losses with yeast orthologs
#
# Jason Jiang - Created: 2022/08/15
#               Last edited: 2022/08/16
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Use k-medoids clustering to cluster microsporidia species by their
#       domain losses w.r.t to their yeast orthologs
#
# Thanks to Brandon Murareanu for this RScript layout
#
# Credit to https://dpmartin42.github.io/posts/r/cluster-mixed-types for the
# inspiration to this script
#
# -----------------------------------------------------------------------------

require(tidyverse)
require(cluster)
require(Rtsne)

################################################################################

## Global variables

set.seed(42)  # for consistent results
PERPLEXITY <- 9  # perplexity value for t-SNE plots, adjust to your liking

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
  microsp_clades <- make_clades_hash(read_csv(args[2], show_col_types = F))
  fam_to_clan <- read_rds(args[3])
  out <- args[4]
  
  species_clan_losses <- get_species_lost_clans(aligned_domain_archs,
                                                fam_to_clan)
  
  # create t-SNE plot of the microsporidia species clustered by their lost
  # domains
  cluster_microsp_by_lost_clans(species_clan_losses,
                                microsp_clades,
                                out)
}

################################################################################

## Helper functions

make_clades_hash <- function(microsp_clades) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  clades <- new.env()
  map2_chr(microsp_clades$species,
           microsp_clades$clade,
           function(x, y) {clades[[x]] <- y})
  
  return(clades)
}


get_species_lost_clans <- function(aligned_domain_archs, fam_to_clan) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    fill_missing_clans(
      aligned_domain_archs %>%
        filter(!is.na(lost_doms), !exclude_species, is_microsp) %>%
        select(species, lost_doms) %>%
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


cluster_microsp_by_lost_clans <- function(species_clan_losses,
                                          microsp_clades,
                                          out) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # calculate pairwise distances between each species based on dissimiliarity
  # in lost clans across orthologs
  gower_dist <- daisy(species_clan_losses[,-1], metric = 'gower')
  
  # select number of clusters, k, to use for k-medoids clustering of the species,
  # based on which k maximizes silhouette width of all clusters
  #
  # TODO - vectorize all of this loopy code
  sil_width <- c(NA)
  for (i in 2 : (nrow(species_clan_losses) - 1)) {
    pam_fit <- pam(gower_dist, diss = T, k = i)
    sil_width[i] <- pam_fit$silinfo$avg.width
  }
  
  # get k that maximizes silhouette width for the clusters
  # k = 2 seems to always generate greatest silhoutte width, which we don't want
  # so, we choose the k after 2 which generates the next greatest silhouette width
  max_sil_width <- which.max(sil_width[3 : length(sil_width)]) + 2
  
  # plot results of k-medoids clustering with t-SNE dimensionality reduction
  plot_lost_dom_tsne(pam(gower_dist, diss = TRUE, k = max_sil_width),
                     gower_dist,
                     species_clan_losses,
                     microsp_clades,
                     out)
}


plot_lost_dom_tsne <- function(k_medoids_clusters, gower_dist,
                               species_clan_losses, microsp_clades,
                               out) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # plot microsporidia species using their dissimilarities in lost clan content,
  # using t-SNE to reduce dimensions to 2
  #
  # note to self: iterate from perplexity = 5 to perplexity = 12
  tsne_obj <- Rtsne(gower_dist, is_distance = T, theta = 0,
                    perplexity = PERPLEXITY)
  
  # note to self: refer to old clustering results for ways to label data, and
  # load that in as a hashmap to color points
  tsne_data <- tsne_obj$Y %>%
    data.frame() %>%
    setNames(c('X', 'Y')) %>%
    mutate(cluster = factor(k_medoids_clusters$clustering),
           species = species_clan_losses$species) %>%
    arrange(cluster) %>%
    rowwise() %>%
    mutate(clade = microsp_clades[[species]])
  
  tsne_plot <- ggplot(aes(x = X, y = Y), data = tsne_data) +
    geom_point(aes(color = clade), size = 6) +
    geom_text(mapping = aes(label = cluster)) +
    # uncomment this line to label each point with their species names as well
    # geom_text(mapping = aes(label = species), nudge_y = 5.5, angle = 30) +
    labs(x = 't-SNE 1', y = 't-SNE 2', title = str_c('Perplexity = ', PERPLEXITY)) +
    theme_bw() +
    theme(plot.title = element_text(size = 18, face = 'bold'),
          axis.title = element_text(size = 14, face = 'bold'),
          axis.text = element_text(size = 14, color = 'black'),
          legend.title = element_text(size = 18, face = 'bold'),
          legend.text = element_text(size = 14))
  
  ggsave(plot = tsne_plot, filename = out, dpi = 600, units = 'in',
         width = 10.99, height = 7.67)
}

################################################################################

main()