######################################################################
# Generate an IGOR J-anchor file from an IgBLAST auxiliary file      #
######################################################################

#==============================================================================
# Preamble
#==============================================================================

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Source packages and auxiliary functions
source("scripts/aux.R")

# Configure input paths
input_path <- snakemake@input[[1]]

# Configure output paths
output_path <- snakemake@output[[1]]

cat("done.\n")

#==============================================================================
# Import data
#==============================================================================

cat("\nImporting auxiliary file...")
col_names <- c("NAME", "FRAME", "TYPE", "COORD")
tab <- suppressMessages(read_tsv(input_path, col_names = col_names))
cat("done.\n")
cat("Aux file dimensions:", dim(tab), "\n")

#==============================================================================
# Convert to anchor file
#==============================================================================

# Convert to anchor file
cat("\nConverting to IGOR anchor format...")
anchor <- tab %>% transmute(NAME = NAME, ANCHOR = COORD + 1)
cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

cat("\nOutput dimensions:", dim(anchor))
cat("\nWriting output table to file...")
write_delim(anchor, output_path, delim=";", col_names = FALSE)
cat("done.\n")
