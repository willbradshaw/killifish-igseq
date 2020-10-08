######################################################################
# Generate an IGOR V-anchor file from IMGT-gapped sequences          #
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
# Run script
#==============================================================================

# Import data
cat("\nReading IMGT-gapped V-sequences...")
v_imgt <- readDNAStringSet(input_path)
cat("done.\n")
cat(length(v_imgt), "sequences read.\n")

# Infer anchors
cat("\nInferring anchor positions of ungapped sequences...")
anchors <- v_imgt %>% subseq(1,209) %>% as.character %>% 
  gsub("\\.", "", .) %>% str_count
cat("done.\n")

# Create anchor table
cat("\nCreating anchor table...")
anchor_tab <- tibble(NAME = names(v_imgt), ANCHOR = anchors)
cat("done.\n")

# Write output
cat("\nOutput dimensions:", dim(anchor_tab))
cat("\nWriting output table to file...")
write_delim(anchor_tab, output_path, delim=";", col_names = FALSE)
cat("done.\n")
