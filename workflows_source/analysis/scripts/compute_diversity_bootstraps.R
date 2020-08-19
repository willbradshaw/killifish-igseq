#######################################################################
# Compute grouped bootstrap distributions from a Change-O table      ##
#######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")
if (!exists("nboot")) nboot <- 2000 # Use standard bootstrap # repeats
if (!exists("uniform")) uniform <- TRUE # Downsample counts by default
if (!exists("max_n")) max_n <- NULL # Downsample counts by default
if (max_n %in% c("None", "NULL")) max_n <- NULL


# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(DUPCOUNT = "n", CONSCOUNT = "n",
            .default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Generate bootstraps
write_log("\nComputing bootstrap table...", newline = FALSE)
bootstraps <- compute_diversity_bootstraps(tab, group_within, group_between,
                                           clone_field, nboot, NULL,
                                           min_n, max_n, uniform)
log_done()

# Write output
write_log("Writing bootstraps to file...")
write_tsv(bootstraps, outpath)
log_done()
