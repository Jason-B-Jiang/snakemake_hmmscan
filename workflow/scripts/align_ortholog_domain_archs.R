# -----------------------------------------------------------------------------
#
# Align ortholog domain architectures
#
# Jason Jiang - Created: 2022/07/29
#               Last edited: 2022/08/08
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
CRH_HEADER = c('query_id', 'match_id', 'score', 'boundaries', 'resolved',
               'cond_evalue', 'indp_evalue')

# Parameters for Needleman-Wunsch alignments of domain architectures
NW_PARAMS = list(MATCH = 0, MISMATCH = -3, GAP = -10, GAPCHAR = '*')

# List of abbreviated outgroup species names
OUTGROUPS <- c('C_eleg', 'D_disc', 'D_reri', 'D_mela', 'H_sapi', 'S_pomb')

################################################################################

main <- function() {
  # ---------------------------------------------------------------------------
  # Command line arguments:
  #   $1 = directory with cath-resolve-hits output files
  #
  #   $2 = filepath to dataframe of all pairwise single-copy orthologs
  #
  #   $3 = filepath to tsv of Pfam families and their clans (from Pfam ftp
  #   site)
  #
  #   $4 = file of microsporidia species w/ poorly sequenced genomes to exclude
  #
  #   $5 = file of essential yeast genes, for annotating orthogroups based on
  #   essentiality of the yeast ortholog
  #
  #   $6 = filepath to save dataframe of aligned domain architectures to
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
  
  # load in lists of microsporidia species w/ poor quality genomes (to exclude
  # from analysis eventually), and essential genes in yeast
  excluded_microsp <- readLines(args[4])
  essential_yeast_genes <- readLines(args[5])
  
  # get hashtable of orthogroups and their ortholog domain architectures
  #
  # if this hashtable was previously constructed and stored in the resources
  # folder, load it in and use it instead of building it from scratch
  if (file.exists(str_c(dirname(dirname(domain_arch_dir)),
                        '/resources/orthogroup_domain_archs.rds'))) {
    orthogroup_domain_archs <- 
      readRDS(str_c(dirname(dirname(domain_arch_dir)),
                    '/resources/orthogroup_domain_archs.rds'))
    
  } else {
    orthogroup_domain_archs <- get_orthogroup_domain_arch_hash(domain_arch_dir,
                                                               fam_to_clan)
    
    saveRDS(orthogroup_domain_archs,
            str_c(dirname(dirname(domain_arch_dir)),
                  '/resources/orthogroup_domain_archs.rds'))
  }
  
  # assign domain architectures to ortholog pairs in orthogroups_df, and
  # align all pairs of single-copy orthologs with Needleman-Wunsch alignment
  aligned_domain_archs <- assign_and_align_domain_archs(orthogroups_df,
                                                        orthogroup_domain_archs,
                                                        fam_to_clan,
                                                        essential_yeast_genes,
                                                        excluded_microsp)
  
  # verify domain losses in outgroups
  outgroups_DA <- get_outgroups_domain_arch_hashtable(aligned_domain_archs,
                                                      OUTGROUPS)
  
  aligned_domain_archs <- aligned_domain_archs %>%
    rowwise() %>%
    mutate(dom_losses_in_outgroups =
             ifelse(!is.na(lost_doms) & is_microsp & species != 'R_allo' &
                      !is.na(aligned_ortholog_domain_archs),
                    verify_dom_losses_in_outgroups(aligned_ortholog_domain_archs,
                                                   orthogroup,
                                                   OUTGROUPS,
                                                   outgroups_DA),
                    NA))
  
  # write modified orthogroups_df dataframe to specified filepath
  write_csv(aligned_domain_archs, args[6])
}

################################################################################

## Helper functions

make_pfam_clan_hashtables <- function(pfam_clans) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  file_name = pfam_clans
  out <- str_c(dirname(dirname(file_name)),  # file to save this hashtable to
               '/resources/pfam_fam_to_clan.rds')
  
  if (file.exists(out)) {
    return(readRDS(out))
  }
  
  pfam_clans <- read_tsv(pfam_clans, show_col_types = F)
  
  # create hashtables mapping both pfam families back to their clans
  fam_to_clan <- new.env()
  Map(function(fam, clan) {fam_to_clan[[fam]] <- ifelse(is.na(clan), fam, clan)},
      pfam_clans$Family_ID, pfam_clans$Clan_name)
  
  # write hashtable to resources folder
  saveRDS(fam_to_clan, out)
  
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
                                          fam_to_clan,
                                          essential_yeast_genes,
                                          excluded_microsp) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    orthogroups_df %>%
      select(-median_diff, -species_p_val, -species_p_val_adj) %>%
      mutate(essential = yeast_ortholog %in% essential_yeast_genes,
             exclude_species = species %in% excluded_microsp) %>%
      rowwise() %>%
      mutate(yeast_domain_archs = get_ortholog_domain_arch(orthogroup,
                                                           yeast_ortholog,
                                                           orthogroup_domain_archs),
             species_domain_archs = get_ortholog_domain_arch(orthogroup,
                                                             species_ortholog,
                                                             orthogroup_domain_archs)) %>%
      ungroup() %>%
      separate(yeast_domain_archs, c('yeast_domain_arch',
                                     'yeast_domain_arch_clan',
                                     'yeast_domain_bounds'), sep = ' ::: ') %>%
      separate(species_domain_archs, c('species_domain_arch',
                                       'species_domain_arch_clan',
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
                                             fam_to_clan)) %>%
      ungroup() %>%
      mutate(lost_doms_len = get_lost_doms_len(aligned_ortholog_domain_archs,
                                               yeast_domain_bounds,
                                               species_domain_arch),
             ortholog_short_enough_for_dom_loss =
               species_len <= (yeast_len - 0.85 * lost_doms_len))
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

### Helper functions for verifying domain losses heuristically

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

## Helper functions for verifying that lost domains in Microsporidia are conserved
## in outgroup species

get_outgroups_domain_arch_hashtable <- function(aligned_domain_archs, outgroups) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  outgroups_DA_hash <- new.env()
  for (sp in outgroups) {
    sp_DA_df <- filter(aligned_domain_archs, species == sp,
                       !is.na(aligned_ortholog_domain_archs))
    sp_orthogroup_domain_archs <- new.env()
    
    Map(function(og, domain_arch) {sp_orthogroup_domain_archs[[og]] <- domain_arch},
        sp_DA_df$orthogroup,
        sp_DA_df$aligned_ortholog_domain_archs)
    
    outgroups_DA_hash[[sp]] <- sp_orthogroup_domain_archs
  }
  
  return(outgroups_DA_hash)
}


verify_dom_losses_in_outgroups <- function(aligned_DA, orthogroup, outgroups,
                                           outgroups_DA_hash) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  # positions of lost domains in the ortholog, relative to the yeast ortholog
  lost_dom_idx <- get_lost_dom_idx(aligned_DA)
  dom_loss_verified <- c()
  
  for (i in lost_dom_idx) {
    dom_loss_verified <- c(dom_loss_verified,
                           is_dom_conserved_in_outgroups(i, orthogroup, outgroups,
                                                         outgroups_DA_hash))
  }
  
  return(str_c(dom_loss_verified, collapse = '; '))
}


get_lost_dom_idx <- function(aligned_DA) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  sp_1_DA <- str_split(str_split(aligned_DA, '; ')[[1]][1], ' -> ')[[1]]
  return(which(sp_1_DA == '*'))
}


is_dom_conserved_in_outgroups <- function(i, orthogroup, outgroups, outgroups_DA_hash) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  is_conserved_in_outgroup <- c()
  for (sp in outgroups) {
    outgroup_DA_alignment <- outgroups_DA_hash[[sp]][[orthogroup]]
    
    if (is.null(outgroup_DA_alignment)) {
      is_conserved_in_outgroup <- c(is_conserved_in_outgroup, FALSE)
      
    } else {
      sp_DA <- str_split(str_split(outgroup_DA_alignment, '; ')[[1]][1], ' -> ')[[1]]
      yeast_DA <- str_split(str_split(outgroup_DA_alignment, '; ')[[1]][2], ' -> ')[[1]]
      yeast_dom_idx <- get_yeast_dom_idx(i, yeast_DA)
      
      is_conserved_in_outgroup <- c(is_conserved_in_outgroup,
                                    sp_DA[yeast_dom_idx] != '*')
    }
  }
  
  return(sum(is_conserved_in_outgroup) >= 1)
}


get_yeast_dom_idx <- function(i, yeast_DA) {
  # ---------------------------------------------------------------------------
  # Docstring goes here.
  # ---------------------------------------------------------------------------
  non_gaps <- 0
  j = 1
  
  while (non_gaps < i) {
    if (yeast_DA[j] == '*') {
      j <- j + 1
    } else {
      non_gaps <- non_gaps + 1
      j <- j + 1
    }
  }
  
  return(j - 1)
}

################################################################################

main()