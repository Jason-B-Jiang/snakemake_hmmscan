# -----------------------------------------------------------------------------
#
# Verify domain losses by ortholog lengths
#
# Jason Jiang - Created: 2022/04/29
#               Last edited: 2022/08/08
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Implement heuristic to verify domain losses inferred from Pfam domain
# assignments for orthologs wrt to yeast, by seeing if the orthologs have gotten
# short enough to have lost their domains.
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

require(tidyverse)

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = dataframe of aligned ortholog domain architectures
  #   (from align_ortholog_domain_archs.R)
  #
  #   $2 = filepath to save annotated aligned ortholog dataframe into
  # ---------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  aligned_domain_archs <- read_csv(args[1], show_col_types = F)
  out <- args[2]
  
  # add new column to aligned_domain_archs indicating whether the non-yeast
  # ortholog has gotten short enough to have truly lost the putatively missing
  # domains from its yeast ortholog
  verified_lost_doms <- verify_lost_doms(aligned_domain_archs)
  
  write_csv(verified_lost_doms, str_c(out, 'aligned_domain_archs.csv'), sep = '/')
}

################################################################################

### Helper functions

verify_lost_doms <- function(aligned_domain_archs) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  return(
    aligned_domain_archs %>%
      mutate(lost_doms_len = get_lost_doms_len(aligned_ortholog_domain_archs,
                                               yeast_domain_bounds,
                                               species_domain_arch),
             ortholog_short_enough_for_dom_loss =
               species_len <= (yeast_len - 0.85 * lost_doms_len))
  )
}


get_lost_doms_len <- Vectorize(function(aligned_DA, yeast_domain_bounds,
                                        species_domain_arch) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  # if no domain assignments for this species ortholog, return the length of
  # all domains in the yeast ortholog
  if (is.na(species_domain_arch)) {
    lost_dom_lens <- str_split(yeast_domain_bounds, '; ')[[1]]
    
  } else {
    # otherwise, get length of all lost domains only from the yeast ortholog
    species_DA <- str_split(str_split(aligned_DA, '; ')[[1]][1], ' -> ')[[1]]
    lost_dom_posns <- which(species_DA == '*')
    
    if (length(lost_dom_posns) == 0) {
      # no lost domains in this species ortholog, so return NA
      return(NA)
    }
    
    lost_dom_lens <- str_split(yeast_domain_bounds, '; ')[[1]][lost_dom_posns]
  }
  
  return(sum(sapply(lost_dom_lens, function(x) {get_domain_len(x)})))
}, vectorize.args = c('aligned_DA', 'yeast_domain_bounds', 'species_domain_arch'))


get_domain_len <- Vectorize(function(dom_bounds) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  dom_fragments <- str_split(dom_bounds, ',')[[1]]
  
  sum(
    sapply(dom_fragments,
           function(x) {as.integer(str_split(x, '-')[[1]][2]) - 
               as.integer(str_split(x, '-')[[1]][1]) + 1})
    )
})

################################################################################

main()