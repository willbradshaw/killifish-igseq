#######################################################################
# Compute diversity spectra for each diversity replicate in dataset   #
#######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Read Change-O table of bootstrap abundances
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(ITER = "i", N = "i", 
            .default = col_character())
bootstraps <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(bootstraps))

# Compute per-replicate diversity spectra
write_log("\nComputing alpha-diversities...", newline = FALSE)
q <- q_range(min_q, max_q, step_q)
alpha_div_bs <- compute_alpha_diversity_bootstrap(bootstraps, group_within,
                                                  group_between, 
                                                  clone_field, q$range)
log_done()

write_log("\nComputing gamma-diversities...", newline = FALSE)
gamma_div_bs <- compute_gamma_diversity_bootstrap(bootstraps, group_within,
                                                  clone_field, q$range)
log_done()

write_log("\nComputing beta-diversities...", newline = FALSE)
beta_div_bs <- compute_beta_diversity_bootstrap(gamma_div_bs, alpha_div_bs, group_within)
log_done()

write_log("\nComputing solo diversities...", newline = FALSE)
solo_div_bs <- compute_solo_diversity_bootstrap(bootstraps, group_within,
                                                  group_between,
                                                  clone_field, q$range)
log_done()

# Write output
write_log("Writing diversity spectra to file...")
write_tsv(alpha_div_bs, outpath_alpha)
write_tsv(beta_div_bs, outpath_beta)
write_tsv(gamma_div_bs, outpath_gamma)
write_tsv(solo_div_bs, outpath_solo)
log_done()
