######################################################################
# Count mutation frequency from a Change-O DB                       ##
######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Set defaults
if (!exists("sequence_field")) sequence_field <- "CLONAL_SEQUENCE"
if (!exists("germline_field")) germline_field <- "CLONAL_GERMLINE"

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Compute observed mutations
write_log("\nComputing observed mutation frequencies...", newline = FALSE)
observed <- observedMutations(tab,
                              sequenceColumn=sequence_field,
                              germlineColumn=germline_field,
                              regionDefinition=IMGT_V,
                              nproc=snakemake@threads)
log_done()
write_log("New dimensions:", dim(observed))

# Write output
write_log("\nSaving to output...", newline = FALSE)
writeChangeoDb(observed, outpath)
log_done()
