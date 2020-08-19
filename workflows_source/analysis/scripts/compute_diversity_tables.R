#######################################################################
# Compute grouped alpha/beta/gamma diversities from a Change-O table ##
#######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")
if (!exists("ci")) ci <- 0.95 # Use standard default confidence intervals
if (!exists("nboot")) nboot <- 2000 # Use standard bootstrap # repeats
if (!exists("uniform")) uniform <- TRUE # Downsample counts by default
if (!exists("max_n")) max_n <- NULL
if (copy == "NULL") copy <- NULL
if (max_n %in% c("None", "NULL")) max_n <- NULL


# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(DUPCOUNT = "n", CONSCOUNT = "n",
            .default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Generate diversity table
write_log("\nComputing diversities...", newline = FALSE)
div <- compute_diversity_tables(tab, group_within, group_between,
                                clone_field, nboot, copy, min_n,
                                max_n, min_q, max_q, step_q, ci,
                                uniform, progress = FALSE)
log_done()

# Write output
write_log("Writing diversity spectrum to file...")
write_tsv(div, outpath)
log_done()
