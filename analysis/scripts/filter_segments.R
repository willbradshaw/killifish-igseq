###################################################################
# Filter a Change-O DB by whether its segment calls are ambiguous #
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
filter_V <- snakemake@params[["filter_V"]]
filter_D <- snakemake@params[["filter_D"]]
filter_J <- snakemake@params[["filter_J"]]

cat("done.\n")

#==============================================================================
# Import data
#==============================================================================

cat("\nImporting Change-O table...")
col <- cols(HAS_V = "l", HAS_D = "l", HAS_J = "l",
            V_AMBIG = "l", D_AMBIG = "l", J_AMBIG = "l",
            .default = col_character())
tab <- suppressMessages(read_tsv(input_path, col_types = col))
cat("done.\n")
cat("Input dimensions:", dim(tab), "\n")

#==============================================================================
# Remove ambiguous segment rows
#==============================================================================

for (segment in c("V", "D", "J")){
    if (get(paste0("filter_", segment))){
        cat("\nFiltering rows with ambiguous", segment, "calls...")
        n <- nrow(tab)
        tab <- filter(tab, pull(tab, paste0("HAS_", segment)),
                      !pull(tab, paste0(segment, "_AMBIG")))
        cat("done.\n")
        n_filtered <- nrow(tab)
        cat("\n   Rows discarded:", as.integer(n - n_filtered))
        cat("\n   New dimensions:", n_filtered)
    } else {
        cat("\nNot filtering on", segment, paste0("calls", "."))
    }
}

#==============================================================================
# Save output
#==============================================================================

cat("\nOutput dimensions:", dim(tab), "\n")
cat("\nWriting output table to file...")
write_tsv(tab, output_path)
cat("done.\n")
