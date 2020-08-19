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
tab <- suppressMessages(read_tsv(inpath_tab, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Read SHM model
write_log("\nImporting mutation model...", newline = FALSE)
mod <- readRDS(inpath_mod)
log_done()

# Count expected mutations and append MU_EXPECTED columns to the output
write_log("\nComputing expected mutation frequencies...", newline = FALSE)
expected <- expectedMutations(tab, 
                              sequenceColumn=sequence_field,
                              germlineColumn=germline_field,
                              targetingModel=mod,
                              regionDefinition=IMGT_V,
                              nproc=snakemake@threads)
log_done()
write_log("New dimensions:", dim(expected))

# Write output
write_log("\nSaving to output...", newline = FALSE)
writeChangeoDb(expected, outpath)
log_done()
