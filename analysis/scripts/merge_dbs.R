#################################################
# Combine entries from two or more Change-O DBs #
#################################################

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
input_paths <- snakemake@input

# Configure output paths
output_path <- snakemake@output[[1]]

cat("done.\n")

#==============================================================================
# Handle input
#==============================================================================

# Check there is at least some input
cat("Number of input paths:", length(input_paths),"\n")
if (length(input_paths) == 0){
    cat("Error: No input files given!")
    stop()
}

# Ignore empty files
file_sizes <- sapply(input_paths, function(p) file.info(p)$size)
input_paths <- input_paths[file_sizes != 0]

# If no non-empty files, stop here
if (length(input_paths) == 0){
    cat("All input files empty! Writing empty file...")
    file.create(output_path) # Pass on empty file
    cat("done.\n")

# If only one non-empty file, just copy that
} else if (length(input_paths) == 1){
    cat("Only one non-empty file; copying to output...")
    col <- cols(.default = col_character())
    db <- suppressMessages(read_tsv(input_paths[[1]], col_types = col))
    db <- mutate_at(db, vars(one_of("DIST_NEAREST")),as.numeric)
    write_tsv(db, output_path)
    cat("done.\n")

# Otherwise, read all files, join them and output
} else {
    cat("Multiple non-empty files; concatenating...")
    col <- cols(.default = col_character())
    dbs <- lapply(input_paths, function(f)
                  suppressMessages(read_tsv(f, col_types = col)))
    # Handle case where all DIST_NEAREST are NA:
    dbs <- lapply(dbs, function(f) mutate_at(f, vars(one_of("DIST_NEAREST")),
                                             as.numeric))
    db <- bind_rows(dbs)
    cat("done.\n")
    cat("\nInput DB rows:", sapply(dbs, nrow))
    cat("\nOutput DB rows:", nrow(db))
    cat("\nDifference in row number:", sum(sapply(dbs, nrow)) - nrow(db))
    cat("\nWriting concatenated database to output file...")
    write_tsv(db, output_path)
    cat("done.\n")
}
