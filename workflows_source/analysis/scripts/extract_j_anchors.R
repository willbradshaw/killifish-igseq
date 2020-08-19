######################################################################
# Generate an IGOR J-anchor file from an IgBLAST auxiliary file      #
######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Read Change-O table
col_names <- c("NAME", "FRAME", "TYPE", "COORD")
write_log("\nImporting auxiliary file...", newline = FALSE)
tab <- suppressMessages(read_tsv(inpath, col_names = col_names))
log_done()
write_log("Aux file dimensions:", dim(tab))

# Convert to anchor file
write_log("\nConverting to IGOR anchor format...", newline = FALSE)
anchor <- tab %>% transmute(NAME = NAME, ANCHOR = COORD + 1)
log_done()

# Write output
write_log("\nOutput dimensions:", dim(anchor))
write_log("\nWriting output table to file...", newline = FALSE)
write_delim(anchor, outpath, delim=";", col_names = FALSE)
log_done()
