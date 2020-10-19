######################################################################
# Prepare sequence DB for sequence extraction for IGoR analysis      #
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

# Configure script parameters
collapsed <- assign_default("collapsed", FALSE)
individuals <- snakemake@params[["individuals"]]
group_field <- snakemake@params[["group_field"]]
groups <- snakemake@params[["groups"]]

cat("done.\n")
cat("Individuals:", individuals, "\n")
cat("Group field:", group_field, "\n")
cat("Groups:", groups, "\n")

#==============================================================================
# Import data
#==============================================================================

cat("\nImporting Change-O table...")
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(input_path, col_types = col))
cat("done.\n")
cat("Input dimensions:", dim(tab), "\n")

#==============================================================================
# Execute script
#==============================================================================

# Filter individuals
cat("\nFiltering individuals...")
tab <- tab %>% filter(INDIVIDUAL %in% as.character(individuals))
cat("done.\n")
cat(nrow(tab), "sequence entries remaining.")

# Filter groups
cat("\nFiltering", group_field, "groups...")
tab <- tab %>% filter(!!as.name(group_field) %in% as.character(groups))
cat("done.\n")
cat(nrow(tab), "sequence entries remaining.")

# Remove gap characters from clonal sequence
cat("\nStripping ambiguous/gap characters from sequences...",
          newline = FALSE)
if (collapsed){
    tab <- tab %>% mutate(SEQUENCE_OUT = gsub("[^ACGTN]", "", CLONAL_SEQUENCE))
} else {
    tab <- tab %>% mutate(SEQUENCE_OUT = gsub("[^ACGTN]", "", SEQUENCE_IMGT))
}
cat("done.\n")
cat(nrow(tab), "sequence entries remaining.")

#==============================================================================
# Save output
#==============================================================================

cat("\nOutput dimensions:", dim(tab))
cat("\nWriting output table to file...")
write_tsv(tab, output_path)
cat("done.\n")
