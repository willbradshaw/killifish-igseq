######################################################
# Compute a pairwise RDI matrix from segment counts ##
######################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Read Change-O table
write_log("\nImporting segment counts...", newline = FALSE)
counts <- read.table(inpath, check.names = FALSE)
log_done()
write_log("Input dimensions:", dim(counts))

if (nrow(counts) == 0){
    write_log("No data; writing empty output file...")
    file.create(outpath)
    log_done()
} else {
    # Generate segment count matrix
    write_log("\nComputing RDI...", newline=FALSE)
    rdi <- calcRDI(counts, subsample = TRUE,
                   distMethod = distance, nIter = iterations,
                   constScale = constant_scale,
                   units = ifelse(transform, "lfc", "pct"))
    log_done()

    # Write output
    write_log("Writing RDI matrix to file...")
    write.csv(as.matrix(rdi), file=outpath, quote = FALSE)
    log_done()

}

# Wildcards: set, segments, group
# Params: distance, iterations, constant_scale, transform, aux
