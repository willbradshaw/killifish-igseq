###############################################################################
## AUXILIARY SCRIPT                                                          ##
## Rarefaction of clonal counts                                              ##
###############################################################################

#==============================================================================
# Preamble
#==============================================================================

# Set up logging
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Source formatting (mainly for packages here)
source("scripts/aux_format-plots.R")

# Set parameters
individuals_excluded <- snakemake@params[["individuals_excluded"]]
n_repeats <- snakemake@params[["n_repeats"]]
sample_sizes <- snakemake@params[["sample_sizes"]]
scale <- snakemake@params[["scale"]]

# Specify input paths
ageing_path <- snakemake@input[["ageing"]]
gut_path <- snakemake@input[["gut"]]

#==============================================================================
# Auxiliary functions
#==============================================================================

get_clonetab <- function(tab, scale, clone_field,
                         group_field){
  tab <- tab %>% filter(!is.na(!!as.name(clone_field)))
  if (is.null(scale)){
    clones <- tab %>% select(!!as.name(group_field), !!as.name(clone_field))
  } else {
    clones <- tibble(N = 1:sum(tab[[scale]]))
    clones[[group_field]] <- rep(tab[[group_field]], tab[[scale]])
    clones[[clone_field]] <- rep(tab[[clone_field]], tab[[scale]])
    clones <- select(clones, -N)
  }
  return(clones)
}

sample_clones_single <- function(tab, sample_size, scale,
                                 clone_field, group_field,
                                 replace){
  # Get individual-labelled vector of clone IDs
  clones <- get_clonetab(tab, scale, clone_field, group_field) %>%
    group_by(!!as.name(group_field)) %>%
    filter(n() >= sample_size)
  # Sample
  samples <- clones %>% sample_n(sample_size, replace = replace)
  # Count clones in new sample
  rarefied <- samples %>% 
    group_by(!!as.name(group_field), !!as.name(clone_field)) %>%
    summarise(N = n()) %>% 
    group_by(!!as.name(group_field)) %>%
    mutate(CLNRANK = row_number(desc(N)), CLNFREQ = N/sum(N),
           SINGLE = N == 1, SMALL = N < 5, LARGE = N >= 5)
  rarefied_summ <- summarise(rarefied, N_CLONES = n(), 
                             N_CLONES_SINGLE = sum(SINGLE),
                             N_CLONES_SMALL = sum(SMALL),
                             N_CLONES_LARGE = sum(LARGE)) %>%
    mutate(PC_CLONES_SMALL = N_CLONES_SMALL/N_CLONES)
  return(rarefied_summ)
}

rarefy_clones_single <- function(tab, sample_size, n_repeats,
                                 scale, clone_field,
                                 group_field, replace){
  # Repeatedly downsample and obtain clone counts for a single
  # sample size
  rarefactions <- lapply(1:n_repeats, function(m)
    sample_clones_single(tab, sample_size, scale, 
                         clone_field, group_field, replace) %>%
    mutate(ITER = m)) %>%
    bind_rows %>% melt(id.vars = c("INDIVIDUAL", "ITER"), 
                       variable.name = "METRIC", value.name = "VALUE")
  rarefactions_summ <- rarefactions %>% 
    group_by(!!as.name(group_field), METRIC) %>%
    summarise(MEAN = mean(VALUE), SD = sd(VALUE)) %>%
    mutate(SAMPLE_SIZE = sample_size)
  return(rarefactions_summ)
}

rarefy_clones <- function(tab, sample_sizes, n_repeats, scale,
                          clone_field = "CLONE",
                          group_field = "INDIVIDUAL",
                          replace = FALSE){
  # Repeatedly downsample and obtain clone counts for a range
  # of sample sizes
  bind_rows(lapply(sample_sizes, function(s)
    rarefy_clones_single(tab, s, n_repeats, scale, clone_field, group_field,
                         replace)))
}

cat("...done.\n")

#==============================================================================
# Import data
#==============================================================================

cat("\nImporting data...")

tab_gut <- import_div(gut_path, 
                      col_types = cols(REPLICATE = "c", FISH = "c",
                                       INDIVIDUAL = "c")) %>%
  filter(!INDIVIDUAL %in% individuals_excluded)

tab_age <- import_div(ageing_path) %>%
  filter(!INDIVIDUAL %in% individuals_excluded)

cat("...done.\n")

#==============================================================================
# Rarefy clonal counts
#==============================================================================

cat("\nRarefying ageing dataset...")
r_age <- rarefy_clones(tab_age, sample_sizes, n_repeats, scale) %>%
  mutate(EXPERIMENT = "ageing")
cat("...done.\n")

cat("\nRarefying gut dataset...")
r_gut <- rarefy_clones(tab_gut, sample_sizes, n_repeats, scale) %>%
  mutate(EXPERIMENT = "gut")
cat("...done.\n")

r <- bind_rows(r_gut, r_age)

#------------------------------------------------------------------------------
# SAVE OUTPUT
#------------------------------------------------------------------------------

cat("\nWriting output...")
write_tsv(r, snakemake@output[[1]])
cat("...done.\n")