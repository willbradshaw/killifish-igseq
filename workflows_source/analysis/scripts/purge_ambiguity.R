#######################################################################
# Convert ambiguous base calls to N                                   #
#######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Specify defaults
if (!exists("sequence_field")) sequence_field <- "SEQUENCE_IMGT"
if (!exists("germline_field")) germline_field <- "GERMLINE_IMGT_D_MASK"

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Remove residual ambiguous characters
chars_allowed <- c("A","T","G","C","N",".","-")
regex_disallowed <- paste0("[^", paste(chars_allowed, collapse = ""), "]")
tab[[sequence_field]] <- gsub(regex_disallowed, "N", tab[[sequence_field]])
tab[[germline_field]] <- gsub(regex_disallowed, "N", tab[[germline_field]])

# Write to output
write_log("\nOutput dimensions:", dim(tab))
write_log("Writing output to file...")
write_tsv(tab, outpath)
log_done()
