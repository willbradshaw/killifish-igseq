################################################################
# Add VJ and VDJ combined segment call fields to a Change-O DB #
################################################################

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
    col <- cols(HAS_V = "l", HAS_D = "l", HAS_J = "l",
                V_AMBIG = "l", D_AMBIG = "l", J_AMBIG = "l",
                .default = col_character())
    tab <- suppressMessages(read_tsv(inpath, col_types = col))
    log_done()
    write_log("Input dimensions:", dim(tab))

    # Copy solo identities from "best" calls
    write_log("\nExtracting V/D/J calls:")
    for (segment in c("V", "D", "J")){
        call <- paste0(segment, "_CALL")
        write_log(call)
        best <- paste0("BEST_", call)
        write_log("   ", "Copying from:", best)
        tab[[call]] <- tab[[best]]
        write_log("   ", "done.")
        write_log("   ", "Output dimensions:", dim(tab))
    }

    # Assign combined identities
    write_log("\nAdding combined segment calls...", newline = FALSE)
    tab <- tab %>% unite(VJ_CALL, V_CALL, J_CALL, sep="/", remove=FALSE) %>%
        unite(VDJ_CALL, V_CALL, D_CALL, J_CALL, sep="/", remove=FALSE) %>%
        mutate(HAS_VJ = HAS_V & HAS_J, HAS_VDJ = HAS_V & HAS_D & HAS_J,
               VJ_AMBIG = V_AMBIG | J_AMBIG,
               VDJ_AMBIG = V_AMBIG | D_AMBIG | J_AMBIG)
    log_done()
    write_log("Output dimensions:", dim(tab))

    # Write output
    write_log("Writing output table to file...")
    writeChangeoDb(tab, outpath)
    log_done()
}
