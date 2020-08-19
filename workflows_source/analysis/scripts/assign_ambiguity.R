######################################################################
# Label rows of a Change-O DB with ambiguity status of segment calls #
######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Test for empty input
if (file.info(inpath)$size == 0){
    write_log("No data; writing empty output file...")
    file.create(outpath)
    log_done()
} else {
    # Read Change-O table
    write_log("\nImporting Change-O table...", newline = FALSE)
    tab <- suppressMessages(readChangeoDb(inpath))
    log_done()
    write_log("Input dimensions:", dim(tab))

    # Generate target and new field names

    write_log("\nGenerating field names...", newline = FALSE)
    prefixes <- c("V","D","J")
    log_done()
    write_log("\nInput fields:", paste0(prefixes, "_CALL"))
    write_log("\nOutput fields:", paste0(prefixes, "_AMBIG"))

    write_log("\nAssigning ambiguities:")
    for (p in prefixes){
        # Identify ambiguous calls and add to tibble
        incol <- paste0(p, "_CALL")
        outcol <- paste0(p, "_AMBIG")
        write_log("   ", incol, "->", outcol)
        ambig <- (is.na(pull(tab, incol)) | str_count(pull(tab, incol), ",") > 0)
        ambig_total <- sum(ambig)/length(ambig)
        write_log("   ", "   ", "% ambiguous:", signif(ambig_total*100))
        tab[[outcol]] <- ambig
    }
    tab <- tab %>% mutate(VJ_AMBIG = V_AMBIG | J_AMBIG,
                          VDJ_AMBIG = V_AMBIG | D_AMBIG | J_AMBIG)

    # Write output
    write_log("\nOutput dimensions:", dim(tab))
    write_log("\nWriting output table to file...")
    writeChangeoDb(tab, outpath)
    log_done()
}
