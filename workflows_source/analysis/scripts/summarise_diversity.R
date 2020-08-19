#########################################################
# Summarise diversity spectra over bootstrap replicates #
#########################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

if (!exists("ci")) ci <- 0.95 # Use standard default confidence intervals
if (divtype != "solo") group_between <- NA
q <- q_range(min_q, max_q, step_q)

# Read Change-O table of bootstrapped diversities
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(ITER = "i", N_GROUP = "i", Q = "d", D = "d",
            .default = col_character())
diversities_bootstrapped <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(diversities_bootstrapped))

# Summarise diversities
write_log("\nComputing summarised diversities...", newline = FALSE)
diversities_summarised <- summarise_bootstrap_diversity(diversities_bootstrapped,
                                                        group_within, ci, q$keep_0,
                                                        group_between) %>%
    mutate(DIVTYPE = divtype)
log_done()

# Write output
write_log("Writing summarised spectra to file...")
write_tsv(diversities_summarised, outpath)
log_done()
