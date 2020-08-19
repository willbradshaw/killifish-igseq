#########################################################
# Summarise diversity spectra over bootstrap replicates #
#########################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Read Change-O table of bootstrapped diversities
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(ITER = "i", N_GROUP = "i", Q = "i", D = "i",
            .default = col_character())
diversities_bootstrapped <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(diversities_bootstrapped))

# Fix groups for solo comparisons
if (divtype == "solo"){
    diversities_bootstrapped[[paste0(group_within, "_OLD")]] <- diversities_bootstrapped[[group_within]]
    diversities_bootstrapped[[group_within]] <- paste(diversities_bootstrapped[[group_within]],
                                                      diversities_bootstrapped[[group_between]],
                                                      sep = "_")
}

# Summarise diversities
write_log("\nComputing pairwise p-values...", newline = FALSE)
diversities_tested <- compute_diversity_pvalues(diversities_bootstrapped, group_within) %>%
    mutate(DIVTYPE = divtype)
log_done()

# Write output
write_log("Writing p-values to file...")
write_tsv(diversities_tested, outpath)
log_done()
