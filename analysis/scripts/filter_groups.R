###################################################################
# Filter a Change-O DB by group specifications                    #
###################################################################

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

# Configure script parameters
filter_groups <- snakemake@params[["filter_groups"]]
keep_complement <- snakemake@params[["keep_complement"]]

cat("done.\n")

#==============================================================================
# Import data
#==============================================================================

cat("\nImporting Change-O table...")
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(snakemake@input[[1]], col_types = col))
cat("done.\n")
cat("Input dimensions:", dim(tab), "\n")

#==============================================================================
# Process data
#==============================================================================

# Filter columns
cat("\nFiltering on specified columns...")
for (col in names(filter_groups)){
    cat("\n   ", col)
    vals_all <- unique(tab[[col]])
    vals_listed <- filter_groups[[col]]
    keep <- vals_all %in% vals_listed
    if (keep_complement) keep <- !keep
    vals_keep <- vals_all[keep]
    tab <- tab[tab[[col]] %in% vals_keep,]
}
cat("\ndone.\n")

#==============================================================================
# Save output
#==============================================================================

# Write output
cat("\nOutput dimensions:", dim(tab), "\n")
cat("\nWriting output table to file...")
write_tsv(tab, output_path)
cat("done.\n")
