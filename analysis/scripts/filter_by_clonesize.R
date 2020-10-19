###################################################################
# Filter a Change-O DB by clone size                              #
###################################################################
# NB: Currently only groups by individual; needs refinement if
# want to group by replicate instead / as well

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
if ("min_size" %in% names(snakemake@params)){
    min_size <- snakemake@params[["min_size"]]
} else {
    min_size <- 1
}
if ("max_size" %in% names(snakemake@params)){
    max_size <- snakemake@params[["max_size"]]
} else {
    max_size <- Inf
}

cat("done.\n")

#==============================================================================
# Import data
#==============================================================================

cat("\nImporting Change-O table...")
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(input_path, col_types = col))
cat("done.\n")
cat("Input dimensions:", dim(tab), "\n")

#==============================================================================
# Count unique sequences per clone per individual
#==============================================================================

cat("\nCounting clone sizes...")
clone_sizes <- tab %>% group_by(INDIVIDUAL, CLONE, SEQUENCE_IMGT) %>%
    summarise(CLNCOUNT = 1) %>% group_by(INDIVIDUAL, CLONE) %>%
    summarise(CLNCOUNT = sum(CLNCOUNT))
tab_cs <- full_join(tab, clone_sizes, by = c("INDIVIDUAL", "CLONE"))
cat("done.\n")
cat("Dimensions:", dim(tab_cs), "\n")

#==============================================================================
# Filter sequences with missing clonal identities
#==============================================================================

cat("\nDiscarding sequences with missing clonal identities...")
tab_nona <- tab_cs %>% filter(!is.na(CLONE))
cat("done.\n")
cat("Dimensions:", dim(tab_nona), "\n")
cat("Entries discarded:", nrow(tab_cs) - nrow(tab_nona), "\n")

#==============================================================================
# Filter sequences by clone size
#==============================================================================

cat("\nDiscarding sequences by clone size...")
tab_sized <- tab_nona %>% filter(between(CLNCOUNT, min_size, max_size))
cat("done.\n")
cat("Dimensions:", dim(tab_sized), "\n")
cat("Entries discarded:", nrow(tab_nona) - nrow(tab_sized), "\n")

#==============================================================================
# Save output
#==============================================================================

cat("\nWriting output table to file...")
write_tsv(tab_sized, output_path)
cat("done.\n")
