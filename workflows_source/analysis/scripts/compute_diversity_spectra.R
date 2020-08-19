#######################################################################
# Compute diversity spectra from pre-computed bootstraps             ##
#######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

if (!exists("ci")) ci <- 0.95 # Use standard default confidence intervals

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(ITER = "i", N = "i", 
            .default = col_character())
bootstraps <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(bootstraps))

# Generate bootstraps
write_log("\nComputing diversity spectra...", newline = FALSE)
div <- compute_diversity_spectra(bootstraps, group_within, group_between,
                                 clone_field, min_q, max_q, step_q, ci)
log_done()

# Write output
write_log("Writing diversity spectra to file...")
write_tsv(div, outpath)
log_done()
