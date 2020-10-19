#################################################
# Combine entries from two or more Change-O DBs #
#################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# ----------------------------------------------------------------------------
# Handle input
# ----------------------------------------------------------------------------

# Handle input paths (independently of aux function)
# TODO: Handle named input files correctly!
write_log("Parsing input paths...", newline = FALSE)
inpaths <- unlist(snakemake@input)
log_done()
write_log(length(inpaths), "input paths given.")
write_log(inpaths)

if (length(inpaths) == 0){
    write_log("Error: No input files given!")
    stop()
}

# Ignore empty files
write_log("\nFiltering empty files...", newline = FALSE)
inpaths <- inpaths[file.info(inpaths)$size != 0]
log_done()
write_log(length(inpaths), "input paths remaining.")

# If no non-empty files, stop here
if (length(inpaths) == 0){
    write_log("All input files empty! Writing empty file...", newline=FALSE)
    file.create(outpath) # Pass on empty file
    log_done()

# If only one non-empty file, just copy that
} else if (length(inpaths) == 1){
    write_log("Only one non-empty file; copying to output...", newline=FALSE)
    file.copy(inpaths[1], outpath, overwrite = TRUE)
    log_done()

# Otherwise, read all files, join them and output
} else {
    col <- cols(.default = col_character())
    write_log("\nMultiple non-empty files; concatenating...", newline=FALSE)
    dbs <- lapply(inpaths, function(f)
                  suppressMessages(read_tsv(f, col_types = col)))
    # Handle case where all of a column is NA by converting all to character
    dbs <- lapply(dbs, function(t) mutate_all(t, as.character))
    db <- bind_rows(dbs) # TODO: Add error handling for nonmatching columns
    log_done()
    write_log("Input DB rows:", sapply(dbs, nrow))
    write_log("Output DB rows:", nrow(db))
    write_log("Difference in row number:", sum(sapply(dbs, nrow)) - nrow(db))
    write_log("\nWriting concatenated database to output file...", newline=FALSE)
    writeChangeoDb(db, outpath)
    log_done()
}
