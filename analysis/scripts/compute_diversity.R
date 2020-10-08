#######################################################################
# Compute diversity spectra for each diversity replicate in dataset   #
#######################################################################

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
output_path_alpha <- snakemake@output[["alpha"]]
output_path_beta <- snakemake@output[["beta"]]
output_path_gamma <- snakemake@output[["gamma"]]
output_path_solo <- snakemake@output[["solo"]]

# Configure script parameters
min_q <- snakemake@params[["min_q"]]
max_q <- snakemake@params[["max_q"]]
step_q <- snakemake@params[["step_q"]]
clone_field <- snakemake@params[["clone_field"]]
group_within <- snakemake@params[["group_within"]]
group_between <- snakemake@params[["group_between"]]

cat("done.\n")

#==============================================================================
# Auxiliary functions
#==============================================================================

compute_solo_diversity_bootstrap <- function(bootstraps, group_within,
                                             group_between, clone_field, Q){
    # Compute solo diversity spectrum within each inner group and bootstrap rep
    bs_tab <- bootstraps %>%
        dplyr::group_by(!!sym(group_within), !!sym(group_between), ITER)
    div_tab_split <- tibble()
    for (q in Q){
        dt <- dplyr::summarise(bs_tab, Q = q, D = calcDiversity(N, q = q),
                        N_GROUP = 1)
        div_tab_split <- bind_rows(div_tab_split, dt)
    }
    return(div_tab_split)
}

compute_gamma_diversity_bootstrap <- function(bootstraps, group_within,
                                              clone_field, Q){
    # Compute gamma diversity within each outer group and bootstrap rep
    bs_tab <- bootstraps %>% 
      dplyr::group_by(!!sym(group_within), !!sym(clone_field), ITER) %>%
      dplyr::summarise(N = sum(N)) %>%
      dplyr::group_by(!!sym(group_within),  ITER)
    div_tab <- tibble()
    for (q in Q){
        dt <- dplyr::summarise(bs_tab, Q = q, D = calcDiversity(N, q = q))
        div_tab <- bind_rows(div_tab, dt)
    }
    return(div_tab)
}

compute_alpha_diversity_bootstrap <- function(bootstraps, group_within,
                                              group_between, clone_field, Q){
    # Compute alpha diversity within each outer group and bootstrap rep
    bs_tab <- bootstraps %>%
      dplyr::group_by(!!sym(group_within), !!sym(group_between), ITER)
    div_tab_split <- tibble()
    for (q in Q){
        dt <- dplyr::summarise(bs_tab, Q = q, D = calcDiversity(N, q = q))
        div_tab_split <- bind_rows(div_tab_split, dt)
    }
    div_tab <- div_tab_split %>% 
        dplyr::group_by(!!sym(group_within), ITER, Q) %>% 
        dplyr::summarise(D = ifelse(dplyr::first(Q) != 1,
                             mean(D^(1-dplyr::first(Q)))^(1/(1-dplyr::first(Q))),
                             exp(mean(log(D)))),
                  N_GROUP = n())
    return(div_tab)
}

compute_beta_diversity_bootstrap <- function(gamma_div_tab, alpha_div_tab,
                                             group_within){
    # Compute beta diversity spectrum within each outer group and bootstrap rep,
    # from pre-computed gamma and alpha spectra
    full_join(gamma_div_tab, alpha_div_tab, by = c(group_within, "ITER", "Q"),
              suffix = c("_GAMMA", "_ALPHA")) %>%
        dplyr::mutate(D= D_GAMMA/D_ALPHA) %>% dplyr::select(-D_GAMMA, -D_ALPHA)
}

#==============================================================================
# Import data
#==============================================================================

cat("\nImporting bootstrap table...")
col <- cols(ITER = "i", N = "i", 
            .default = col_character())
bootstraps <- suppressMessages(read_tsv(input_path, col_types = col))
cat("done.\n")
cat("Input dimensions:", dim(bootstraps), "\n")

#==============================================================================
# Compute diversities from bootstraps
#==============================================================================

# Compute per-replicate diversity spectra
cat("\nComputing alpha-diversities...")
q <- q_range(min_q, max_q, step_q)
alpha_div_bs <- compute_alpha_diversity_bootstrap(bootstraps, group_within,
                                                  group_between, 
                                                  clone_field, q$range)
cat("done.")

cat("\nComputing gamma-diversities...")
gamma_div_bs <- compute_gamma_diversity_bootstrap(bootstraps, group_within,
                                                  clone_field, q$range)
cat("done.")

cat("\nComputing beta-diversities...")
beta_div_bs <- compute_beta_diversity_bootstrap(gamma_div_bs, alpha_div_bs, group_within)
cat("done.")

cat("\nComputing solo diversities...")
solo_div_bs <- compute_solo_diversity_bootstrap(bootstraps, group_within,
                                                  group_between,
                                                  clone_field, q$range)
cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

cat("\nWriting diversity spectra to file...")
write_tsv(alpha_div_bs, output_path_alpha)
write_tsv(beta_div_bs, output_path_beta)
write_tsv(gamma_div_bs, output_path_gamma)
write_tsv(solo_div_bs, output_path_solo)
cat("done.\n")
