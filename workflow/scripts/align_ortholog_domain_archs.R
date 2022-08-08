# -----------------------------------------------------------------------------
#
# Align ortholog domain architectures
#
# Jason Jiang - Created: 2022/07/29
#               Last edited: 2022/08/07
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

require(tidyverse)
require(NameNeedle)

################################################################################

## Define global variables

# Column names for cath-resolve-hits output files
CRH_HEADER = c('query_id', 'match_id', 'score', 'boundaries',
               'resolved', 'cond_evalue', 'indp_evalue')

# Parameters for Needleman-Wunsch alignments of domain architectures
NW_PARAMS = list(MATCH = 0, MISMATCH = -3, GAP = -10, GAPCHAR = '*')

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = directory with cath-resolve-hits output files
  #   $2 = filepath to dataframe of all pairwise single-copy orthologs
  #   $3 = filepath to tsv of Pfam families and their clans (from Pfam ftp
  #        site)
  #   $4 = filepath to save dataframe of aligned domain architectures to
  # ---------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = T)
  
  # set directory to read in cath-resolve-hits domain architecture files from
  # for each orthogroup
  domain_arch_dir <- args[1]
  
  # load in dataframe of orthogroups with single copy yeast gene and at least
  # single copy microsporidia ortholog, from Orthofinder
  orthogroups_df <- read_csv(args[2], show_col_types = F)
  
  # make hashtable mapping Pfam domain families back to their clans
  fam_to_clan <- make_pfam_clan_hashtables(args[3])
  
  # get hashtable of orthogroups and their ortholog domain architectures
  orthogroup_domain_archs <- get_orthogroup_domain_arch_hash(domain_arch_dir,
                                                             fam_to_clan)
  
  # assign domain architectures to ortholog pairs in orthogroups_df, and
  # align all pairs of single-copy orthologs with Needleman-Wunsch alignment
  aligned_domain_archs <- assign_and_align_domain_archs(orthogroups_df,
                                                        orthogroup_domain_archs,
                                                        fam_to_clan)
  
  # write modified orthogroups_df dataframe to specified filepath
  write_csv(aligned_domain_archs, args[4])
}

################################################################################

## Helper functions

make_pfam_clan_hashtables <- function(pfam_clans) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  file_name = pfam_clans
  pfam_clans <- read_tsv(pfam_clans, show_col_types = F)
  
  # create hashtables mapping both pfam families back to their clans
  fam_to_clan <- new.env()
  Map(function(fam, clan) {fam_to_clan[[fam]] <- ifelse(is.na(clan), fam, clan)},
      pfam_clans$Family_ID, pfam_clans$Clan_name)
  
  # write hashtable to resources folder
  saveRDS(fam_to_clan, str_c(dirname(dirname(file_name)),
                             '/resources/pfam_fam_to_clan.rds'))
  
  return(fam_to_clan)
}


get_orthogroup_domain_arch_hash <- function(domain_arch_dir, fam_to_clan) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # only consider cath-resolve-hits files that are non-empty
  files = list.files(domain_arch_dir, full.names = T)[
    file.size(list.files(domain_arch_dir, full.names = T)) > 0
  ]
  
  orthogroups = unname(sapply(files, function(x) {basename(x)}))

  domain_arch_hash <- new.env()
  Map(function(og, resolved_hits) {
      domain_arch_hash[[og]] <- get_og_domain_archs(resolved_hits, fam_to_clan)
      },
    orthogroups,
    files)
  
  return(domain_arch_hash)
}


get_og_domain_archs <- function(resolved_hits, fam_to_clan) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # get dataframe of resolved domain architectures from cath-resolve-hits
  # using try-catch block to see cases where we get a warning, just doing
  # as a sanity check for myself
  tryCatch(
    expr = {crh_df <- parse_crh_output(resolved_hits, fam_to_clan)},
    warning = function(w) {
      message(resolved_hits)
    })
  
  # create hashtable mapping each ortholog in the orthogroup to a named list of
  # its domain architecture and domain boundaries
  og_domain_archs <- new.env()
  sapply(1 : nrow(crh_df),
         function(i) {
           og_domain_archs[[crh_df$query_id[i]]] <- list(
             'domain_arch' = crh_df$match_id[i],
             'domain_arch_clan' = crh_df$match_clans[i],
             'domain_bounds' = crh_df$resolved[i]
           )
         })
  
  return(og_domain_archs)
}


parse_crh_output <- function(resolved_hits, fam_to_clan) {
  # ---------------------------------------------------------------------------
  # Read in cath-resolve-hits output file for an orthogroup's ortholog domain
  # architectures + domain boundaries
  #
  # Args:
  #   resolved_hits: filepath to cath-resolve-hits output file
  #
  # ---------------------------------------------------------------------------
  return(
    read_delim(resolved_hits, delim = ' ', comment = '#',
                       col_names = CRH_HEADER, show_col_types = F) %>%
    select(query_id, match_id, resolved) %>%
    group_by(query_id) %>%
    mutate(match_id = str_c(match_id, collapse = '; '),
           resolved = str_c(resolved, collapse = '; ')) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(match_clans = get_domain_clans(match_id, fam_to_clan)) %>%
    ungroup() %>%
    group_by(query_id) %>%
    distinct(.keep_all = T)
  )
}


get_domain_clans <- function(domains, fam_to_clan) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  domains <- str_split(domains, '; ')[[1]]
  
  return(str_c(
    # if domain isn't in fam_to_clan hashtable, just keep the domain as is
    sapply(domains, function(x) {ifelse(!is.null(fam_to_clan[[x]]),
                                        fam_to_clan[[x]],
                                        x)}),
    collapse = '; '))
}

################################################################################

## Helper functions for aligning domain architectures and getting domain
## architecture differences between ortholog pairs

assign_and_align_domain_archs <- function(orthogroups_df,
                                          orthogroup_domain_archs,
                                          fam_to_clan) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    orthogroups_df %>%
      select(-median_diff, -species_p_val, -species_p_val_adj) %>%
      rowwise() %>%
      mutate(yeast_domain_archs = get_ortholog_domain_arch(orthogroup,
                                                           yeast_ortholog,
                                                           orthogroup_domain_archs),
             species_domain_archs = get_ortholog_domain_arch(orthogroup,
                                                             species_ortholog,
                                                             orthogroup_domain_archs)) %>%
      ungroup() %>%
      separate(yeast_domain_archs, c('yeast_domain_arch', 'yeast_domain_arch_clan',
                                     'yeast_domain_bounds'), sep = ' ::: ') %>%
      separate(species_domain_archs, c('species_domain_arch', 'species_domain_arch_clan',
                                       'species_domain_bounds'), sep = ' ::: ') %>%
      rowwise() %>%
      mutate(aligned_ortholog_domain_archs = align_domain_archs(species_domain_arch,
                                                                species_domain_arch_clan,
                                                                yeast_domain_arch,
                                                                yeast_domain_arch_clan),
             lost_doms = get_lost_doms(aligned_ortholog_domain_archs,
                                       species_domain_arch,
                                       yeast_domain_arch),
             gained_doms = get_gained_doms(aligned_ortholog_domain_archs,
                                           species_domain_arch,
                                           yeast_domain_arch),
             swapped_doms = get_swapped_doms(aligned_ortholog_domain_archs,
                                             fam_to_clan))
  )
}


get_ortholog_domain_arch <- function(orthogroup, ortholog, orthogroup_domain_archs) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  if (is_null(orthogroup_domain_archs[[orthogroup]]) | 
      is_null(orthogroup_domain_archs[[orthogroup]][[ortholog]])) {
    return(NA)
  }
  
  return(str_c(
    orthogroup_domain_archs[[orthogroup]][[ortholog]]$domain_arch,
    orthogroup_domain_archs[[orthogroup]][[ortholog]]$domain_arch_clan,
    orthogroup_domain_archs[[orthogroup]][[ortholog]]$domain_bounds,
    sep = ' ::: '
  ))
}


align_domain_archs <- function(DA_1, DA_1_clans, DA_2, DA_2_clans) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  if (any(is.na(DA_1), is.na(DA_2))) {
    return(NA)
  }
  
  DA_1 <- str_split(DA_1, '; ')[[1]]
  DA_1_clans <- str_split(DA_1_clans, '; ')[[1]]
  DA_2 <- str_split(DA_2, '; ')[[1]]
  DA_2_clans <- str_split(DA_2_clans, '; ')[[1]]
  
  # since we are aligning domain architectures using their domain clans, get
  # a list mapping each clan to a unique letter representation
  clans_to_letters <- get_clan_letters(unique(c(DA_1_clans, DA_2_clans)))
  
  # create string representations of each domain architecture using the clan to
  # letter mappings, for alignment
  DA_1_str = str_c(sapply(DA_1_clans, function(x) {clans_to_letters[[x]]}),
                   collapse = '')

  DA_2_str = str_c(sapply(DA_2_clans, function(x) {clans_to_letters[[x]]}),
                   collapse = '')
  
  # align domain architectures using their string representations
  aligned_DA <- needles(DA_1_str, DA_2_str, params = NW_PARAMS)
  
  # return formatted string of domain architecture alignments, with original
  # domains replacing the letter representations in the alignments
  return(
    str_c(
      format_alignment_str(aligned_DA$align1, DA_1),
      format_alignment_str(aligned_DA$align2, DA_2),
      sep = '; '
    )
  )
}


get_clan_letters <- function(clans) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  clans_to_letters <- lapply(as.list(1 : length(clans)), function(i) {LETTERS[i]})
  names(clans_to_letters) <- clans
  
  return(clans_to_letters)
}


format_alignment_str <- function(aln_str, domain_arch) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  aln_str <- str_split(aln_str, '')[[1]]
  
  i = 1
  for (j in 1 : length(aln_str)) {
    if (aln_str[j] != '*') {
      aln_str[j] <- domain_arch[i]
      i <- i + 1
    }
  }
  
  return(str_c(aln_str, collapse = ' -> '))
}


get_lost_doms <- function(aligned_DA, DA_1, DA_2) {
  # ---------------------------------------------------------------------------
  # Return domains from the first ortholog domain architecture, DA_1, that have
  # been lost from the second ortholog domain architecture, DA_2, using their
  # aligned domain architectures.
  #
  # Args:
  #   aligned_DA: aligned domain architectures for both orthologs
  #   
  #   DA_1: domain architecture for first ortholog
  #
  #   DA_2: domain architecture for second ortholog.
  #
  # ---------------------------------------------------------------------------
  # no domain assignments for first ortholog, so just treat all domains from the
  # second ortholog as 'lost' domains
  if (is.na(DA_1)) {
    return(DA_2)
    
    # no domain assignments for second ortholog, so no domain losses possible for
    # ortholog
  } else if (is.na(DA_2)) {
    return(NA)
    
  }
  
  aligned_DA <- str_split(aligned_DA, '; ')[[1]]
  aligned_DA_1 <- str_split(aligned_DA[1], ' -> ')[[1]]
  aligned_DA_2 <- str_split(aligned_DA[2], ' -> ')[[1]]
    
  lost_doms <- aligned_DA_2[which(aligned_DA_1 == '*')]
    
  return(ifelse(length(lost_doms) > 0,
                str_c(lost_doms, collapse = '; '),
                NA))
}


get_gained_doms <- function(aligned_DA, DA_1, DA_2) {
  # ---------------------------------------------------------------------------
  # Return domains from the first ortholog domain architecture, DA_1, that have
  # been gained from the second ortholog domain architecture, DA_2, using their
  # aligned domain architectures.
  #
  # Args:
  #   aligned_DA: aligned domain architectures for both orthologs
  #   
  #   DA_1: domain architecture for first ortholog
  #
  #   DA_2: domain architecture for second ortholog
  #
  # ---------------------------------------------------------------------------
  # if no domain assignments for first ortholog, then no domain gains possible
  # for first ortholog
  if (is.na(DA_1)) {
    return(NA)
    
    # if no domain assignments for second ortholog, then consider all domains in
    # first ortholog to be gained domains
  } else if (is.na(DA_2)) {
    return(DA_1)
    
  }
  
  aligned_DA <- str_split(aligned_DA, '; ')[[1]]
  aligned_DA_1 <- str_split(aligned_DA[1], ' -> ')[[1]]
  aligned_DA_2 <- str_split(aligned_DA[2], ' -> ')[[1]]
    
  gained_doms <- aligned_DA_1[which(aligned_DA_2 == '*')]
    
  return(ifelse(length(gained_doms) > 0,
                str_c(gained_doms, collapse = '; '),
                NA))
}


get_swapped_doms <- function(aligned_DA, fam_to_clan) {
  # ---------------------------------------------------------------------------
  # Return domains from the first ortholog domain architecture, DA_1, that have
  # been swapped from the second ortholog domain architecture, DA_2, using their
  # aligned domain architectures.
  #
  # Args:
  #   aligned_DA: aligned domain architectures for both orthologs
  #   
  #   DA_1: domain architecture for first ortholog
  #
  #   DA_2: domain architecture for second ortholog
  #
  # ---------------------------------------------------------------------------
  # if no domain assignments for either ortholog, then no domain swapping is
  # possible between the two orthologs
  if (is.na(aligned_DA)) {
    return(NA)
  }
  
  aligned_DA <- str_split(aligned_DA, '; ')[[1]]
  
  # convert domains in domain architectures back to their clans
  # domain swap = domain of different CLAN in ortholog 1 than in ortholog 2
  aligned_DA_1 <- str_split(aligned_DA[1], ' -> ')[[1]]
  aligned_DA_1_clans <- sapply(aligned_DA_1,
                               function(x) {ifelse(x != '*',
                                                   fam_to_clan[[x]],
                                                   x)})
  
  aligned_DA_2 <- str_split(aligned_DA[2], ' -> ')[[1]]
  aligned_DA_2_clans <- sapply(aligned_DA_2,
                              function(x) {ifelse(x != '*',
                                                  fam_to_clan[[x]],
                                                  x)})
  
  swapped_doms_idx <- which(aligned_DA_1_clans != aligned_DA_2_clans & 
                              aligned_DA_1_clans != '*' &
                              aligned_DA_2_clans != '*')
  
  swapped_doms <- map2_chr(aligned_DA_1[swapped_doms_idx],
                           aligned_DA_2[swapped_doms_idx],
                           function(x, y) {str_c(x, y, sep = ' -> ')})
  
  return(ifelse(length(swapped_doms) > 0,
                str_c(swapped_doms, collapse = '; '),
                NA))
}

################################################################################

main()