# -----------------------------------------------------------------------------
#
# Plot distribution of relative lengths of orthologs to yeast genes
#
# Jason Jiang - Created: 2022/08/22
#               Last edited: 2022/08/22
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: For each species, plot the distribution of gene lengths with respect
# to the lengths of their yeast orthologs.
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# ------------------------------------------------------------------------------

require(tidyverse)
require(svglite)

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = csv with single-copy orthologs from all species to yeast, and their
  #   lengths
  #
  #   $2 = file with list of all essential yeast genes, to annotate ortholog
  #   pairs as being yeast-essential or not
  #
  #   $3 = directory to save resulting plots to
  # ---------------------------------------------------------------------------
  orthologs <- read_csv(args[1], show_col_types = F)
  essential_yeast_genes <- read_lines(args[2])
  out_dir <- args[3]
  
  # create the output directory
  dir.create(out_dir)
  
  # annotate whether each species-yeast ortholog pair is essential or not
  orthologs <- orthologs %>%
    mutate(essential = yeast_ortholog %in% essential_yeast_genes)
  
  # stratify dataframe by species
  microsp_orthologs <- filter(orthologs, is_microsp, species != 'R_allo')
  rozella_orthologs <- filter(orthologs, species == 'R_allo')
  outgroup_orthologs <- filter(orthologs, !is_microsp)
  
  # create plots for each individual group of species, and save
  Map(function(orthologs, name, out) {plot_ortholog_length_distn(orthologs, name, out_dir)},
      list(microsp_orthologs, rozella_orthologs, outgroup_orthologs),
      list('All microsporidia', 'Rozella', 'Outgroup species'))
}

################################################################################

## Helper functions

get_length_bin <- Vectorize(function(length_ratio) {
  # ---------------------------------------------------------------------------
  # Get the bin that the ratio of a species ortholog's length to its yeast
  # ortholog falls into.
  #
  # Arguments:
  #   length_ratio: length of ortholog of species divided by length of its
  #   yeast ortholog
  # ---------------------------------------------------------------------------
  if (length_ratio < 0.3) {
    return('<30%')
  } else if (0.3 <= length_ratio & length_ratio <= 0.6) {
    return('30% - 60%')
  } else if (0.6 <= length_ratio & length_ratio <= 0.9) {
    return('60% - 90%')
  } else {
    return('>90%')
  }
})


plot_ortholog_length_distn <- function(orthologs, name, out) {
  # ---------------------------------------------------------------------------
  # Plot the distribution of relative lengths of a species' genes to its
  # yeast single-copy orthologs.
  #
  # Arguments:
  #   orthologs: dataframe of ortholog lengths
  #
  #   name: the group of species the orthologs dataframe represents
  #
  #   out: directory to save resulting plot to
  # ---------------------------------------------------------------------------
  orthologs <- orthologs %>%
    mutate(length_bin = get_length_bin(species_len / yeast_len))
  
  # turn length_bin column into ordered factors, for plotting purposes
  orthologs$length_bin <- factor(orthologs$length_bin,
                                 levels = c('<30%', '30% - 60%', '60% - 90%', '>90%'))
  
  # make and save plot of distribution of ortholog lengths
  length_distn_plot <- ggplot(orthologs, aes(x = length_bin, fill = essential)) +
    geom_bar(position = 'dodge', color = 'black') +
    # add labels for counts on top of bars
    geom_text(aes(label = ..count..), stat = 'count',
              position = position_dodge(0.9), vjust=-0.2) +
    scale_fill_manual(values = c('white', 'black')) +
    labs(x = 'Relative length of ortholog to yeast',
         title = str_c(name, '\nn = ', nrow(orthologs), ' single-copy orthologs to yeast')) +
    theme(title = element_text(size = 18, color = 'black'),
          axis.title = element_text(size = 18, color = 'black'),
          axis.text = element_text(size = 14, color = 'black')) +
    theme_bw()
  
  ggsave(str_c(out, '/', name, '.svg'),
         plot=length_distn_plot,
         units='in',
         dpi=600,
         width=5.34,
         height=4.71)
}

################################################################################

main()