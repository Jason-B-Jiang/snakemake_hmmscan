# -----------------------------------------------------------------------------
#
# Plot domain architecture conservation rates between yeast and orthologs from
# Microsporidia + other species
#
# Jason Jiang - Created: 2022/01/25
#               Last edited: 2022/08/08
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: For each species, plot the frequencies of domain loss, gain and swap 
# in their single-copy orthologs to yeast
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# ------------------------------------------------------------------------------

require(tidyverse)
require(svglite)

################################################################################

main <- function() {
  # ----------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = filepath to aligned ortholog domain architectures, produced by
  #   align_ortholog_domain_archs.R
  #
  #   $2 = file of microsporidia species w/ poorly sequenced genomes to exclude
  #
  #   $3 = file of essential yeast genes, for annotating orthogroups based on
  #   essentiality of the yeast ortholog
  #
  #   $4 = output directory to save plots in
  # ----------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  aligned_domain_archs <- read_csv(args[1], show_col_types = F)
  excluded_microsp <- readLines(args[2])
  essential_yeast_genes <- readLines(args[3])
  out_dir <- args[4]
  
  # format aligned domain archs dataframe so microsporidia species with low quality
  # genomes are excluded from the dataframe, all ortholog pairs are annotated with
  # the essentiality of their yeast orthologs and 
  aligned_domain_archs <- format_aligned_domain_archs_df(aligned_domain_archs,
                                                         essential_yeast_genes,
                                                         excluded_microsp)
  
  # split dataframe by species and plot domain arch conservation rates individually
  # for each species
  dir.create(out_dir)  # initialize the directory to write plots into
  sp_domain_archs <- split(aligned_domain_archs, f=aligned_domain_archs$species)
  Map(function(sp_df, sp_name) {draw_DA_conservation_plot(sp_df, sp_name,
                                                          out_dir)},
      sp_domain_archs,
      names(sp_domain_archs))
  
  # draw plot of overall domain arch conservation rates for Microsporidia only
  # TODO - maybe draw plots for each clade of microsporidia?
  draw_DA_conservation_plot(filter(aligned_domain_archs, is_microsp,
                                   species != 'R_allo'),
                            'All microsporidia', out_dir)
  
  # draw plot of overall domain arch conservation rates for non-Microsporidia
  draw_DA_conservation_plot(filter(aligned_domain_archs, !is_microsp),
                            'Outgroup species', out_dir)
}

################################################################################

# Helper functions

format_aligned_domain_archs_df <- function(aligned_domain_archs,
                                           essential_yeast_genes,
                                           excluded_microsp) {
  # ----------------------------------------------------------------------------
  # Docstring goes here.
  # ----------------------------------------------------------------------------
  return(
    aligned_domain_archs %>%
      filter(!(species %in% excluded_microsp)) %>%
      # add column indicating if yeast ortholog in ortholog pair is essential
      mutate(essential = yeast_ortholog %in% essential_yeast_genes,
             # add columns indicating whether a particular microsporidia ortholog
             # has undergone domain gain, swap and/or loss
             # use these columns for domain arch conservation rates later
             lost = as.integer(!is.na(lost_doms)),
             gain = as.integer(!is.na(gained_doms)),
             swap = as.integer(!is.na(swapped_doms)),
             DA_conservation = get_DA_conservation(lost, gain, swap)) %>%
      separate_rows(DA_conservation, sep = ',') %>%
      mutate(essential = ifelse(essential, 'Essential', 'Non-essential'))
    )
}

get_DA_conservation <- Vectorize(function(lost, gain, swap) {
  # ----------------------------------------------------------------------------
  # Docstring goes here.
  # ----------------------------------------------------------------------------
  changes <- c()
  
  if (lost > 0) {
    changes <- append(changes, 'Loss')
  }
  
  if (gain > 0) {
    changes <- append(changes, 'Gain')
  }
  
  if (swap > 0) {
    changes <- append(changes, 'Swap')
  }
  
  if (length(changes) == 0) {
    return('Conserved')
  }
  
  return(str_c(changes, collapse = ','))
}, vectorize.args = c('lost', 'gain', 'swap'))


draw_DA_conservation_plot <- function(sp_df, sp_name, out_dir) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  # format sp_df for plotting, by summarising counts of domain arch change events
  # across essential and non-essential ortholog pairs
  num_ortholog_pairs <- nrow(sp_df)  # get this number before we transform the dataframe
  
  sp_df <- sp_df %>%
    group_by(DA_conservation, essential) %>%
    summarise(n = n())
  
  # plot and save in out_dir directory
  plot <- ggplot(data = sp_df, aes(x = DA_conservation, y = n, fill = essential)) +
    geom_bar(stat = 'identity') +
    labs(x = 'Domain architecture similarity to yeast',
         y = 'count',
         title = str_c(sp_name, '\n',
                       str_c('n = ', num_ortholog_pairs, ' single-copy ortholog pairs'))) +
    theme_bw() +
    theme(legend.position = c(0.73, 0.85),
          axis.title = element_text(size = 12, color = 'black'),
          axis.text = element_text(size = 12, color = 'black'),
          legend.text = element_text(size = 12),
          legend.title = element_blank())
  
  ggsave(str_c(out_dir, '/', sp_name, '.svg'),
         plot=plot,
         units='in',
         dpi=600,
         width=3.83,
         height=3.91)
}

################################################################################

main()  # run the script