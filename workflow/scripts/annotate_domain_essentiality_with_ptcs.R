# -----------------------------------------------------------------------------
#
# Annotate essentiality/dispensablility of lost ortholog domains using yeast
# PTC data
#
# Jason Jiang - Created: 2022/09/19
#               Last edited: 2022/09/20
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

library(tidyverse)
library(readxl)

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # args <- commandArgs(trailingOnly = T)
  sgd_to_uniprot_names <- read_rds('../../data/yeast_genome_to_uniprot.rds')
  ptc_data <- format_ptc_data(readxl::read_xls('../../data/supp_11.xls'),
                              sgd_to_uniprot_names)
  yeast_domain_archs <- format_yeast_domain_archs(
    read_delim('../../resources/yeast_domains_resolved.txt')
  )
  
  # add domain arch data to ptc data
  ptc_data <- left_join(ptc_data, yeast_domain_archs, by = 'gene')
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
             protein_length = protein_length - 1,
             # replace D/A classifications of PTCs w/ survival = FALSE and survival =
             # TRUE
             ptc_tolerance = ifelse(!is.na(ptc_tolerance),
                                    ifelse(ptc_tolerance == 'A', T, F),
                                    NA)) %>%
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
  ptc_tolerance <- as.logical(str_split(ptc_tolerance, ', ')[[1]])
  ptc_position <- str_split(ptc_position, ', ')[[1]]
  
  # corresponding ptc position in gene for last lethal ptc
  critical_ptc_position <- ptc_position[tail(which(!ptc_tolerance), n = 1)]
  
  if (length(critical_ptc_position) == 0) {
    # all PTCs in this gene are tolerated
    return(NA)
  }
  
  return(critical_ptc_position)
}, vectorize.args = c('ptc_tolerance', 'ptc_position'))


format_yeast_domain_archs <- function(yeast_domain_archs) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    yeast_domain_archs %>%
      group_by(`query-id`) %>%
      mutate(`match-id` = str_c(`match-id`, collapse = ', '),
             resolved = str_c(resolved, collapse = '; ')) %>%
      distinct(`query-id`, .keep_all = T) %>%
      ungroup() %>%
      mutate(domain_starts = get_domain_starts_or_ends(resolved, FALSE),
             domain_ends = get_domain_starts_or_ends(resolved, TRUE)) %>%
      # select(`query-id`, `match-id`, domain_starts, domain_ends) %>%
      rename(gene = `query-id`, domain_arch = `match-id`) %>%
      select(gene, domain_arch, domain_starts, domain_ends)
  )
}


get_domain_starts_or_ends <- Vectorize(function(domain_bounds, get_ends) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  domain_bounds <- str_split(domain_bounds, '; ')[[1]]
  
  if (get_ends) {
    # get stop boundaries of each domain in the protein
    bounds <- unname(sapply(
      domain_bounds,
      function(x) {str_split(tail(str_split(x, ',')[[1]], n = 1), '-')[[1]][2]}
    ))
    
  } else {
    # get start boundaries of each domain in the protein
    bounds <- unname(sapply(
      domain_bounds,
      function(x) {str_split(head(str_split(x, ',')[[1]], n = 1), '-')[[1]][1]}
    ))
  }
  
  return(str_c(bounds, collapse = ', '))
}, vectorize.args = c('domain_bounds', 'get_ends'))


annotate_domains_with_ptc_data <- function(ptc_tolerance, ptc_position,
                                           domain_starts, domain_ends,
                                           critical_ptc_position) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  if (is.na(domain_starts) | is.na(domain_ends)) {
    # NA values for domain starts/ends = no domain assignments for this gene
    return(NA)
  }
  
  # convert inputs into right types
  ptc_tolerance <- as.logical(str_split(ptc_tolerance, ', '))
  ptc_position <- as.integer(str_split(ptc_position, ', '))
  domain_starts <- as.integer(str_split(domain_starts, ', '))
  domain_ends <- as.integer(str_split(domain_ends, ', '))
  critical_ptc_position <- as.integer(critical_ptc_position)
  
  return(
    str_c(
      map2_chr(domain_starts, domain_ends,
               function(domain_start, domain_end) {
                 classify_domain(domain_start,
                                 domain_end,
                                 ptc_tolerance,
                                 ptc_position,
                                 critical_ptc_position)
               }),
      collapse = ', '
    )
  )
}


classify_domain <- function(domain_start, domain_end, ptc_tolerance,
                            ptc_position, critical_ptc_position) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  if (any(ptc_position < domain_start & ptc_tolerance)) {
    # Non-lethal PTCs come before the start of the domain, so dispensable
    return('Dispensable')
    
  } else {
    # index of first lethal PTC disrupting this domain
    first_disrupting_ptc <- tail(which(ptc_position < domain_end & !ptc_tolerance),
                                 n = 1)
    
    # Essential domain:
    # 1) lethal PTC comes before domain end, and all PTCs after the first lethal
    #    PTC disrupting the domain are tolerated
    #
    # OR
    #
    # 2) lethal PTC coming before the domain end is the last PTC in the gene
    #
    # OR
    #
    # 3) lethal PTC coming before the domain end is the second to last PTC in
    #    the gene
    if ((first_disrupting_ptc < length(ptc_position) & 
        all(ptc_tolerance[first_disrupting_ptc : length(ptc_tolerance)])) |
        (length(ptc_position) == first_disrupting_ptc | length(ptc_position) == first_disrupting_ptc + 1)) {
          return('Essential')
    }
    
    return('Ambiguous')
  }
}

################################################################################

main()