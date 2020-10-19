#####################################################
# Disambiguate clone labels from different sources ##
#####################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")
if (!exists("clone_field")) clone_field <- "CLONE"
if (!exists("source_field")) source_field <- "SAMPLE"

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log(nrow(tab), "sequence entries imported.")

# Relabel clones
write_log("\nRelabelling clones...", newline = FALSE)
tab[[clone_field]] <- paste(tab[[source_field]],
                            tab[[clone_field]], sep="_")
log_done()

# Write output
write_log("Writing table to file...")
write_tsv(tab, outpath)
log_done()
