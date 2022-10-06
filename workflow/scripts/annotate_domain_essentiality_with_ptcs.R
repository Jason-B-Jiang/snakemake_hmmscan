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

# Ordering of microsporidia/outgroup species
SP_ORDER <- c('A_alge', 'A_locu', 'E_aedi', 'E_bien', 'E_brev', 'E_canc',
              'E_cuni', 'E_hell', 'E_hepa', 'E_inte', 'E_roma', 'H_erio',
              'N_ausu', 'N_cera', 'N_cide', 'N_disp', 'N_ferr', 'N_gran',
              'N_homo', 'N_iron', 'N_majo', 'N_mino', 'N_okwa', 'N_pari',
              'O_coll', 'P_epip', 'P_phil', 'S_loph', 'T_cont', 'T_homi',
              'T_rati', 'V_corn', 'V_culi', 'A_sp.', 'M_incu', 'M_daph',
              'P_sacc', 'R_allo', 'C_eleg', 'D_disc', 'D_mela', 'D_reri',
              'H_sapi', 'S_pomb')

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
  
  # add domain arch data to ptc data, and add column indicating essentiality
  # of each domain in yeast
  ptc_data <- left_join(ptc_data, yeast_domain_archs, by = 'gene') %>%
    rowwise() %>%
    mutate(domain_ptc_annotations = annotate_domains_with_ptc_data(ptc_tolerance,
                                                                   ptc_position,
                                                                   domain_starts,
                                                                   domain_ends,
                                                                   critical_ptc_position))
  
  # get aligned ortholog yeast domain architectures to annotate lost domains
  aligned_domain_archs <- read_csv('../../results/aligned_domain_archs.csv')
  
  lost_doms <- aligned_domain_archs %>%
    filter(!is.na(lost_doms)) %>%
    rename(gene = yeast_ortholog) %>%
    left_join(ptc_data, by = 'gene') %>%
    rename(yeast_ortholog = gene) %>%
    filter(!is.na(domain_ptc_annotations)) %>%
    mutate(lost_dom_annotations = get_lost_dom_annotations(aligned_ortholog_domain_archs,
                                                           domain_ptc_annotations)) %>%
    separate_rows(lost_dom_annotations, sep = ', ')
  
  # note to self: show absolute COUNTS of lost domain types instead per species
  # show that despite Microsporidia sharing less (?) orthogroups on average to
  # yeast than outgroup species, they still lose more of every domain
  plot_lost_dom_types(lost_doms, aligned_domain_archs)
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
  ptc_tolerance <- as.logical(str_split(ptc_tolerance, ', ')[[1]])
  ptc_position <- as.integer(str_split(ptc_position, ', ')[[1]])
  domain_starts <- as.integer(str_split(domain_starts, ', ')[[1]])
  domain_ends <- as.integer(str_split(domain_ends, ', ')[[1]])
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
  if (any(ptc_position <= domain_start & ptc_tolerance)) {
    # Non-lethal PTCs come before the start of the domain, so dispensable
    return('Dispensable')
    
    # All PTCs tolerated in this gene, but none of the PTCs come before the
    # start of this gene - ambiguous
  } else if (is.na(critical_ptc_position)) {
    return('Ambiguous')
    
    # Essential:
    # 1) First lethal PTC in the gene comes before the domain end
    # OR
    # 2) First lethal PTC in the gene comes right after the domain end
  } else if (critical_ptc_position <= domain_end |
             (length(ptc_position[ptc_position > domain_end]) == 1 &&
              ptc_position[ptc_position > domain_end] == critical_ptc_position)) {
    return('Essential')
  }
  
  return('Ambiguous')
}


get_lost_dom_annotations <- Vectorize(function(aligned_ortholog_domain_archs,
                                               domain_ptc_annotations) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ortholog lost all domains in the yeast ortholog
  if (is.na(aligned_ortholog_domain_archs)) {
    return(domain_ptc_annotations)
  }
  
  domain_ptc_annotations <- str_split(domain_ptc_annotations, ', ')[[1]]
  
  aligned_ortholog_domain_archs <- lapply(
    as.list(str_split(aligned_ortholog_domain_archs, '; ')[[1]]),
    function(s) {str_split(s, ' -> ')[[1]]}
  )
  
  return(
    str_c(domain_ptc_annotations[
      which(aligned_ortholog_domain_archs[[1]] == '*')
      ],
        collapse = ', ')
  )
}, vectorize.args = c('aligned_ortholog_domain_archs', 'domain_ptc_annotations'))


plot_lost_dom_types <- function(lost_doms, aligned_domain_archs) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  ggplot(filter(lost_doms, ortholog_short_enough_for_dom_loss), aes(x = species)) +
    geom_bar(aes(fill = lost_dom_annotations)) +
    coord_flip() +
    theme_bw() +
    labs(y = 'Number of lost domains in orthologs to yeast') +
    scale_fill_manual('Lost domain classification',
                      values = c('black', '#008080', '#DC143C')) +
    scale_x_discrete(limits = rev(SP_ORDER)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = 'black', size = 14),
          axis.title.x = element_text(color = 'black', size = 18),
          axis.text.x = element_text(color = 'black', size = 14),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18))
  
  filter(aligned_domain_archs, !is.na(species_domain_arch) | !is.na(yeast_domain_arch)) %>%
    group_by(species) %>%
    summarise(n = n()) %>%
    ggplot(aes(x = species, y = n)) +
    geom_bar(stat = 'identity', fill = 'black') +
    coord_flip() +
    theme_bw() +
    labs(y = 'Single-copy orthologs to yeast\n(with domain assignments)') +
    scale_x_discrete(limits = rev(SP_ORDER)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = 'black', size = 14),
          axis.title.x = element_text(color = 'black', size = 14),
          axis.text.x = element_text(color = 'black', size = 14))
}

################################################################################

# main()