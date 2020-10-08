#########################################################
# Summarise diversity spectra over bootstrap replicates #
#########################################################

#==============================================================================
# Preamble
#==============================================================================

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Source packages and auxiliary functions
source("scripts/aux.R")

# Configure input paths
input_path <- snakemake@input[[1]]

# Configure output paths
output_path <- snakemake@output[[1]]

# Configure script parameters
min_q <- snakemake@params[["min_q"]]
max_q <- snakemake@params[["max_q"]]
step_q <- snakemake@params[["step_q"]]
ci <- snakemake@params[["ci"]]
group_within <- snakemake@params[["group_within"]]
group_between <- snakemake@params[["group_between"]]
divtype <- snakemake@wildcards[["divtype"]]

cat("done.\n")

#==============================================================================
# Auxiliary functions
#==============================================================================

summarise_bootstrap_diversity <- function(div_tab, group_within, ci, keep_0,
                                          group_between = NA,
                                          empirical_bounds = TRUE){
  #' dplyr::summarise diversity table across bootstrap repeats
  # Calculate empirical quantiles
  quantile_upper <- 1-(1-ci)/2
  quantile_lower <- (1-ci)/2
  # Compute summary table
  summ_tab <- div_tab %>% dplyr::group_by(!!sym(group_within), Q) %>%
    {if ("N_GROUP" %in% colnames(.)) dplyr::group_by(., N_GROUP,
                                              add = TRUE) else . } %>%
    {if (!is.na(group_between)) dplyr::group_by(., !!sym(group_between),
                                         add = TRUE) else . } %>%
    dplyr::summarise(D_SD = sd(D),
              D_LOWER_EMPIRICAL = quantile(D, quantile_lower),
              D_UPPER_EMPIRICAL = quantile(D, quantile_upper),
              D = mean(D)
              ) %>%
    dplyr::mutate(D_LOWER_NORMAL = D - ci_scale(ci) * D_SD,
           D_UPPER_NORMAL = D + ci_scale(ci) * D_SD)
    # Add evenness information
    n = ifelse(!is.na(group_between), 2, 1)
    summ_tab_evenness <- summ_tab %>%
        dplyr::group_by(!!sym(group_within)) %>%
        {if (!is.na(group_between)) dplyr::group_by(., !!sym(group_between),
                                              add = TRUE) else . } %>%
        filter(Q==0) %>% transmute(D_0 = D) %>%
        full_join(summ_tab,
                  by = c(group_within, group_between)[1:n]) %>%
        dplyr::mutate(E = D/D_0, 
               E_LOWER_EMPIRICAL = D_LOWER_EMPIRICAL/D_0, 
               E_UPPER_EMPIRICAL = D_UPPER_EMPIRICAL/D_0,
               E_LOWER_NORMAL = D_LOWER_NORMAL/D_0, 
               E_UPPER_NORMAL = D_UPPER_NORMAL/D_0
               ) %>%
        dplyr::select(-D_0)
    # Drop 0-order row if not in desired Q-range
    if (keep_0) return(summ_tab_evenness)
    return(sum_tab_evenness %>% filter(Q != 0))
}

ci_scale <- function(ci){
  # Compute z-score and SD scale from confidence interval, assuming gaussian
  ci_z <- ci + (1 - ci)/2
  return(qnorm(ci_z))
}

#==============================================================================
# Import data
#==============================================================================

cat("\nImporting bootstrapped diversities...")
col <- cols(ITER = "i", N_GROUP = "i", Q = "d", D = "d",
            .default = col_character())
diversities_bootstrapped <- suppressMessages(read_tsv(input_path, col_types = col))
cat("done.\n")
cat("Input dimensions:", dim(diversities_bootstrapped), "\n")

#==============================================================================
# Summarise diversities
#==============================================================================

if (divtype != "solo") group_between <- NA
q <- q_range(min_q, max_q, step_q)

# Summarise diversities
cat("\nComputing summarised diversities...")
diversities_summarised <- summarise_bootstrap_diversity(diversities_bootstrapped,
                                                        group_within, ci, q$keep_0,
                                                        group_between) %>%
    mutate(DIVTYPE = divtype)
cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

cat("Writing summarised spectra to file...")
write_tsv(diversities_summarised, output_path)
cat("done.\n")
