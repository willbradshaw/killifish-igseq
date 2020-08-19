#############################################################
# Compute segment counts across groups for RDI computation ##
#############################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log(nrow(tab), "sequence entries imported.")

if (nrow(tab) == 0){
    write_log("No data; writing empty output file...")
    file.create(outpath)
    log_done()
} else {
    # Generate segment count matrix
    write_log("\nGenerating segment count matrix...", newline=FALSE)
    call_name <- paste0(toupper(segments), "_CALL")
    counts <- calcVDJcounts(genes = pull(tab, call_name),
                            seqAnnot = as.character(pull(tab, group)),
                            simplifyNames = TRUE,
                            splitCommas = FALSE)
    log_done()

    # Write output
    write_log("Writing count table to file...")
    write.table(counts, file=outpath, quote = FALSE)
    log_done()

}
