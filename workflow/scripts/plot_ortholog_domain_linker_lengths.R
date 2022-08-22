# -----------------------------------------------------------------------------
#
# Plot ortholog domain and linker lengths
#
# Jason Jiang - Created: 2022/08/22
#               Last edited: 2022/08/22
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Create boxplots comparing differences in linker and domain lengths 
# between orthologs that species share with yeast.
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# -----------------------------------------------------------------------------

require(tidyverse)
require(ggsignif)

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = csv of species + yeast domain architectures, with domain lengths
  #   $2 = directory to save resulting plots to
  # ---------------------------------------------------------------------------
  aligned_domain_archs <- read_csv(args[1], show_col_types = F)
  out_dir <- args[2]
  
  # create output directory if not already existing
  dir.create(out_dir)
  
  # stratify domain architectures by species
  microsp_domain_archs <- filter(aligned_domain_archs, is_microsp, species != 'R_allo')
  rozella_domain_archs <- filter(aligned_domain_archs, species == 'R_allo')
  outgroup_domain_archs <- filter(aligned_domain_archs, !is_microsp)
  
  # create boxplots comparing linker and domain lengths between ortholog pairs
  # for each group of species w/ yeast
  Map(function(domain_archs, name, out) {
    plot_domain_and_linker_lengths(domain_archs, name, out)},
      list(microsp_domain_archs, rozella_domain_archs, outgroup_domain_archs),
      list('All microsporidia', 'Rosella', 'Outgroup species'))
}

################################################################################

## Helper functions

plot_domain_and_linker_lengths <- function(domain_archs, name, out) {
  # ---------------------------------------------------------------------------
  # Create boxplots comparing linker vs domain lengths between a species'
  # orthologs to yeast.
  #
  # Arguments:
  #   domain_archs:
  #   name:
  #   out:
  #
  # ---------------------------------------------------------------------------
  domain_archs <- get_orthololog_domain_and_linker_lengths(domain_archs)
  
  # p-values for difference between species + yeast ortholog linker and domain
  # lengths
}


get_ortholog_domain_and_linker_lengths <- function(domain_archs) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    domain_archs %>%
      filter(!is.na(aligned_ortholog_domain_archs)) %>%
      select(species, yeast_ortholog, species_ortholog,
             yeast_len, species_len,
             yeast_domain_bounds, species_domain_bounds,
             aligned_ortholog_domain_archs) %>%
      rowwise() %>%
      mutate(yeast_domain_len = get_domain_lengths(yeast_domain_bounds),
             species_domain_len = get_domain_lengths(species_domain_bounds)) %>%
      ungroup() %>%
      mutate(yeast_linker_len = get_linker_length(yeast_len, yeast_domain_len),
             species_linker_len = get_linker_length(species_len, species_domain_len),
             aligned_yeast_domain_len = get_aligned_domain_lengths(aligned_ortholog_domain_archs,
                                                                   yeast_domain_len,
                                                                   T),
             aligned_species_domain_len = get_aligned_domain_lengths(aligned_ortholog_domain_archs,
                                                                     species_domain_len,
                                                                     F)) %>%
      separate_rows(aligned_yeast_domain_len,
                    aligned_species_domain_len,
                    sep = '; ') %>%
      mutate(aligned_yeast_domain_len = as.integer(aligned_yeast_domain_len),
             aligned_species_domain_len = as.integer(aligned_species_domain_len))
  )
}


sum_domain_fragment_lengths <- function(domain_fragments) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    sum(unname(
      sapply(domain_fragments,
             function(d) {
               as.integer(str_split(d, '-')[[1]])[2] - as.integer(str_split(d, '-')[[1]])[1] + 1
             }
      )
    )
    ))
}


get_domain_lengths <- function(domain_len) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    as.list(str_split(domain_len, '; ')[[1]]) %>%
      lapply(function(x) {sum_domain_fragment_lengths(str_split(x, ',')[[1]])}) %>%
      str_c(collapse = '; ')
  )
}


get_linker_length <- Vectorize(function(prt_len, domain_len) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  domain_len <- sum(as.integer(str_split(domain_len, '; ')[[1]]))
  
  # linker length is all residues in the protein not belonging to a domain
  return(prt_len - domain_len)
}, vectorize.args = c('prt_len', 'domain_len'))


get_aligned_domain_lengths <- Vectorize(function(aligned_ortholog_domain_arch,
                                       domain_bounds,
                                       is_yeast) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  i <- ifelse(is_yeast, 2, 1)
  
  aligned_ortholog_domain_arch <- str_split(
    str_split(aligned_ortholog_domain_arch, '; ')[[1]], ' -> '
  )
  
  domain_bounds <- str_split(domain_bounds, '; ')[[1]]
  
  names(aligned_ortholog_domain_arch[[i]]) <-
    integer(length(aligned_ortholog_domain_arch[[i]]))
  
  j = 1
  for (k in 1 : length(aligned_ortholog_domain_arch[[i]])) {
    if (aligned_ortholog_domain_arch[[i]][k] != '*') {
      names(aligned_ortholog_domain_arch[[i]])[k] <- j
      j <- j + 1
    } else {
      names(aligned_ortholog_domain_arch[[i]])[k] <- NA
    }
  }
    
  aligned_posns <- which(aligned_ortholog_domain_arch[[1]] != '*' &
                           aligned_ortholog_domain_arch[[2]] != '*')
  
  return(
    str_c(domain_bounds[as.integer(names(aligned_ortholog_domain_arch[[i]])[aligned_posns])],
          collapse = '; ')
  )
}, vectorize.args = c('aligned_ortholog_domain_arch', 'domain_bounds', 'is_yeast'))

################################################################################

main()