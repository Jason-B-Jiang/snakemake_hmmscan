# -----------------------------------------------------------------------------
#
# Align ortholog domain architectures
#
# Jason Jiang - Created: 2022/07/29
#               Last edited: 2022/08/04
#
# Reinke Lab - Microsporidia Orthologs Project
#
# Goal: Align domain architectures between yeast proteins + orthologs from
# other species
#
# Thanks to Brandon Murareanu for this RScript layout
#
#
# ------------------------------------------------------------------------------

library(tidyverse)

################################################################################

## Define global variables

# Column names for cath-resolve-hits output files
CRH_HEADER = c('query_id', 'match_id', 'score', 'boundaries',
               'resolved', 'cond_evalue', 'indp_evalue')


################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = directory with cath-resolve-hits output files
  #   $2 = filepath to dataframe of all pairwise single-copy orthologs
  #   $3 = filepath to save dataframe of aligned domain architectures to
  # ---------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  
  # set directory to read in cath-resolve-hits domain architecture files from
  # for each orthogroup
  domain_arch_dir <- args[1]
  
  # load in dataframe of orthogroups with single copy yeast gene and at least
  # single copy microsporidia ortholog, from Orthofinder
  orthogroups_df <- read_csv(args[2], show_col_types = F)
  
  # create hashtable mapping each orthogroup in orthogroups_df to the domain
  # architectures of orthologs in each orthogroups
  orthogroup_domain_archs <- make_orthogroup_domain_archs_hash(
    unique(orthogroups_df$orthogroup), domain_arch_dir
    )
  
  # add in columns for domain architectures + domain architecture alignments
  # for each ortholog with yeast
  orthogroups_df <- orthogroups_df %>%
    rowwise() %>%
    mutate(species_domain_arch = get_ortholog_domain_arch(orthogroup,
                                                          species_ortholog,
                                                          orthogroup_domain_archs),
           yeast_domain_arch = get_ortholog_domain_arch(orthogroup,
                                                        yeast_ortholog,
                                                        orthogroup_domain_archs)) %>%
    filter(!is.na(species_domain_arch), !is.na(yeast_domain_arch)) %>%
    ungroup() %>%
    mutate(aligned_domain_archs = align_ortholog_domain_archs(species_ortholog,
                                                              yeast_ortholog))
  
  # write the resulting dataframe of domain architectures + alignments for
  # each orthogroups as a csv
  write_csv(
    orthogroup_df,
    '../../results/resolved_domain_architectures/orthogroup_domain_archs.csv')
}

################################################################################

## Helper functions

make_orthogroup_domain_archs_hash <- function(orthogroups, domain_arch_dir) {
  # ---------------------------------------------------------------------------
  # Create a hashtable mapping orthogroups to hashtables of their orthologs
  # and domain architectures + domain boundaries
  #
  # Ex: {orthogroup: {ortholog_1: c(domain architecture, domain bounds), ...}}
  #
  # Args:
  #   orthogroups_df:
  #
  #   domain_arch_dir:
  #
  # ---------------------------------------------------------------------------
  orthogroup_domain_archs <- new.env()
  
  for (og in orthogroups) {
    og_domain_arch_df <- parse_crh_output(og, domain_arch_dir)
    
    # there are domain assignments for at least one ortholog in this orthogroup
    if (!is.na(og_domain_arch_df)) {
      orthogroup_domain_archs[[og]] <-
        get_orthogroup_domain_archs(og_domain_arch_df)
    }
  }
  
  return(orthogroup_domain_archs)
}


parse_crh_output <- function(og, domain_arch_dir) {
  # ---------------------------------------------------------------------------
  # Read in cath-resolve-hits output file for an orthogroup's (og) ortholog
  # domain architectures + domain boundaries
  #
  # Args:
  #   og:
  #
  #   domain_arch_dir:
  #
  # ---------------------------------------------------------------------------
  og_domain_arch_df <- read_delim(str_c(domain_arch_dir, og, sep = '/'),
                       delim = ' ', comment = '#')
  
  # if no domain assignments at all for this orthogroup, just return NA
  if (nrow(og_domain_arch_df) < 1) {
    return(NA)
  }
  
  colnames(og_domain_arch_df) <- CRH_HEADER
  
  og_domain_arch_df <- og_domain_arch_df %>%
    select(query_id, match_id, resolved) %>%
    group_by(query_id) %>%
    mutate(match_id = str_c(match_id, collapse = '|'),
           resolved = str_c(resolved, collapse = '|')) %>%
    distinct(.keep_all = )
  
  return(og_domain_arch_df)
}


get_orthogroup_domain_archs <- function(og_domain_arch_df) {
  # ---------------------------------------------------------------------------
  # Docstring goes here
  # ---------------------------------------------------------------------------
  # turn this dataframe into a hashtable, with keys as ortholog names and
  # values as vectors of domain architecture strings and domain boundaries
  domain_arch_hash <- new.env()
  
  for (i in 1 : nrow(og_domain_arch_df)) {
    domain_arch_hash[[og_domain_arch_df$query_id[i]]] <-
      c(og_domain_arch_df$match_id[i], og_domain_arch_df$resolved[i])
  }
  
  return(domain_arch_hash)
}


get_ortholog_domain_arch <- function(orthogroup, ortholog,
                                     orthogroup_domain_archs) {
  # ---------------------------------------------------------------------------
  # Docstring goes here
  # ---------------------------------------------------------------------------
  domain_archs <- orthogroup_domain_archs[[orthogroup]][[ortholog]]
  if (is.null(domain_lengths)) {
    return(NA)
  }
  
  return(domain_archs[1])
}


get_ortholog_domain_lengths <- function(orthogroup, ortholog,
                                        orthogroup_domain_archs) {
  # ---------------------------------------------------------------------------
  # Docstring goes here
  # ---------------------------------------------------------------------------
  domain_lengths <- orthogroup_domain_archs[[orthogroup]][[ortholog]]
  if (is.null(domain_lengths)) {
    return(NA)
  }
  
  return(domain_lengths[2])
}


align_ortholog_domain_archs <- function(species_ortholog, yeast_ortholog) {
  
}