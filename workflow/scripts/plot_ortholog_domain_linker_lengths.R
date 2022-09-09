# -----------------------------------------------------------------------------
#
# Plot ortholog domain and linker lengths
#
# Jason Jiang - Created: 2022/08/22
#               Last edited: 2022/09/02
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
  args <- commandArgs(trailingOnly = T)
  aligned_domain_archs <- read_csv(args[1], show_col_types = F) %>%
    filter(!exclude_species)  # remove orthologs from species w/ low qual genomes
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
    plot_domain_and_linker_lengths(domain_archs, name, out_dir)},
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
  domain_archs <- get_ortholog_domain_and_linker_lengths(domain_archs)
  
  # dataframe just for unique linker values
  linkers <- distinct(domain_archs, species_ortholog, .keep_all = T) %>%
    select(species, essential, yeast_ortholog, species_ortholog,
           species_linker_len, yeast_linker_len, linker_ratio) %>%
    rename(length_ratio = linker_ratio,
           species_len = species_linker_len,
           yeast_len = yeast_linker_len) %>%
    mutate(type = 'linker') %>%
    group_by(essential) %>%
    mutate(sd_by_essential = sd(length_ratio, na.rm=T),
           mean_by_essential = mean(length_ratio, na.rm=T),
           is_outlier = abs(length_ratio - mean(domains$length_ratio, na.rm = T)) > 1.5 * sd_by_essential)
  
  # dataframe just for aligned domain lengths
  domains <- domain_archs %>%
    select(species, essential, yeast_ortholog, species_ortholog,
           aligned_species_domain_len, aligned_yeast_domain_len,
           domain_ratio) %>%
    rename(length_ratio = domain_ratio,
           species_len = aligned_species_domain_len,
           yeast_len = aligned_yeast_domain_len) %>%
    mutate(type = 'domain') %>%
    group_by(essential) %>%
    mutate(sd_by_essential = sd(length_ratio, na.rm=T),
           mean_by_essential = mean(length_ratio, na.rm=T),
           is_outlier = abs(length_ratio - mean(domains$length_ratio, na.rm = T)) > 1.5 * sd_by_essential)
  
  # p-values for difference between species + yeast ortholog linker and domain
  # lengths
  p_domains <- wilcox.test(domains$species_len,
                           domains$yeast_len,
                           paired = T)$p.value  # -4 AA med diff, 97% of length
  p_linkers <- wilcox.test(linkers$species_len,
                           linkers$yeast_len,
                           paired = T)$p.value
  
  # p-values for differences in domain/linker lengths between essential and
  # non-essential ortholog pairs
  
  # 97% in essential, 95.5% in non-essential
  # greater significance than linkers due to larger sample size
  p_domains_essential <- wilcox.test(filter(domains, essential)$length_ratio,
                                     filter(domains, !essential)$length_ratio)$p.value
  
  # 73.6% in essential, 71.3% in non-essential
  p_linkers_essential <- wilcox.test(filter(linkers, essential)$length_ratio,
                                     filter(linkers, !essential)$length_ratio)$p.value
  
  # do some really weird stuff for plotting
  dom_linkers <- rbind(domains, linkers)  # combine dataframes for plotting
  
  dom_plot <- ggplot(data = filter(dom_linkers, type == 'domain', !is_outlier),
                     aes(x = essential, y = length_ratio)) +
    geom_violin(show.legend = F, aes(fill = essential)) +
    stat_summary(fun = "median", fun.min = "median", fun.max= "median",
                 size= 0.3, geom = "crossbar") +
    geom_signif(comparisons = list(c('TRUE', 'FALSE'))) +
    labs(x = 'Essential?',
         y = 'ln(domain length in species ortholog / domain length in yeast ortholog)',
         title = str_c(name, '\nn = ', nrow(linkers),
                       ' species-yeast ortholog pairs\n')) +
    theme(axis.title = element_text(color = 'black', size = 24),
          axis.text = element_text(color = 'black', size = 24),
          plot.title = element_text(size = 24)) +
    theme_bw()
  
  linker_plot <- ggplot(data = filter(dom_linkers, type == 'linker', !is_outlier),
                         aes(x = essential, y = length_ratio)) +
    geom_violin(show.legend = F, aes(fill = essential)) +
    stat_summary(fun = "median", fun.min = "median", fun.max= "median",
                 size= 0.3, geom = "crossbar") +
    geom_signif(comparisons = list(c('TRUE', 'FALSE'))) +
    labs(x = 'Essential?',
         y = 'ln(linker length in species ortholog / linker length in yeast ortholog)',
         title = str_c(name, '\nn = ', nrow(linkers),
                       ' species-yeast ortholog pairs\n')) +
    theme(axis.title = element_text(color = 'black', size = 24),
          axis.text = element_text(color = 'black', size = 24),
          plot.title = element_text(size = 24)) +
    theme_bw()
  
  # helper function to round p values to n significant digits
  # taken from https://stackoverflow.com/questions/43050903/round-to-significant-digits-only-with-decimal-portion-of-number-in-r
  my_signif = function(x, digits) floor(x) + signif(x %% 1, digits)
  
  # dataframe holding overall species vs yeast domain and linker length diffs
  # p values
  # f_labels <- data.frame(type = c('domain', 'linker'),
  #                        label = c(str_c('overall: p = ', my_signif(p_domains, 2)),
  #                                  str_c('overall: p = ', my_signif(p_linkers, 2))))
  
  # add overall p-values for microsporidia domains vs yeast domains lengths,
  # microsporidia linkers vs yeast linkers lengths to the boxplot
  # boxplt <- boxplt +
  #   geom_label(x = 1.9, y = -5.5, aes(label = label), data = f_labels)
  
  # save the plots to desired output directory
  ggsave(plot = dom_plot, filename = str_c(out, '/', name, ' domains', '.svg'),
         units = 'in', width = 5.38, height = 6.88, dpi = 600)
  
  ggsave(plot = linker_plot, filename = str_c(out, '/', name, ' linkers', '.svg'),
         units = 'in', width = 5.38, height = 6.88, dpi = 600)
}


get_ortholog_domain_and_linker_lengths <- function(domain_archs) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  return(
    domain_archs %>%
      filter(!is.na(aligned_ortholog_domain_archs)) %>%
      select(species, essential, yeast_ortholog, species_ortholog,
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
             aligned_species_domain_len = as.integer(aligned_species_domain_len),
             domain_ratio = aligned_species_domain_len / aligned_yeast_domain_len,
             linker_ratio = ifelse(yeast_linker_len != 0,
                                   species_linker_len / yeast_linker_len,
                                   NA))
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