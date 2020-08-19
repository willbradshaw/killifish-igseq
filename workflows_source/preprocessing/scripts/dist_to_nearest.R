#################################################################
# Compute nearest-neighbour Hamming distances for a Change-O DB #
#################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")
separate_failed <- exists("outpath_f")
write_log("Saving dropped sequences separately:", separate_failed)

# ----------------------------------------------------------------------------
# Test for empty input by reading first line from file, then import
# ----------------------------------------------------------------------------

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col)) %>%
    mutate(ROW = seq(n()))
log_done()
write_log(nrow(tab), "sequence entries imported.\n")

# ----------------------------------------------------------------------------
# Compute nearest-neighbour distances
# ----------------------------------------------------------------------------

write_log("Computing nearest-neighbour distances...", newline = FALSE)
dist_pass <- distToNearest(tab, model = model, normalize = normalize,
                      symmetry = symmetry, nproc = snakemake@threads,
                      fields = group, first = FALSE)
log_done()
write_log(nrow(dist_pass), "sequence entries in output.")
write_log("Proportion of rows with nearest-neighbour distances:", 
    mean(!is.na(dist_pass$DIST_NEAREST)), "\n")

# Retrieve dropped sequences
write_log("Retrieving dropped sequences...", newline = FALSE)
dist_fail <- tab %>% filter(! ROW %in% dist_pass$ROW) %>% 
    mutate(DIST_NEAREST = NA)
log_done()
write_log(nrow(dist_fail), "dropped sequence entries retrieved.")
write_log(nrow(dist_fail) + nrow(dist_pass), "sequence entries total.")
write_log("Proportion of total (pass+fail) rows with nearest-neighbour distances:",
    mean(!is.na(c(dist_pass$DIST_NEAREST, dist_fail$DIST_NEAREST))), "\n")

# ----------------------------------------------------------------------------
# Save output database
# ----------------------------------------------------------------------------

if (separate_failed){
    # Write output
    write_log("Writing table of passed sequences to file...", newline = FALSE)
    writeChangeoDb(dist_pass, outpath_p)
    log_done()
    write_log("Writing table of dropped sequences to file...", newline = FALSE)
    writeChangeoDb(dist_fail, outpath_f)
    log_done()
} else {
    write_log("Writing table of all sequences to file...", newline = FALSE)
    writeChangeoDb(bind_rows(dist_pass, dist_fail), outpath)
    log_done()
}
