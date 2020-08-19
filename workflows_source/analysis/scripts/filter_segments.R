###################################################################
# Filter a Change-O DB by whether its segment calls are ambiguous #
###################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(HAS_V = "l", HAS_D = "l", HAS_J = "l",
            V_AMBIG = "l", D_AMBIG = "l", J_AMBIG = "l",
            .default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Remove ambiguous segment rows
for (segment in c("V", "D", "J")){
    if (get(paste0("filter_", segment))){
        write_log("\nFiltering rows with ambiguous", segment, "calls...",
                  newline = FALSE)
        n <- nrow(tab)
        tab <- filter(tab, pull(tab, paste0("HAS_", segment)),
                      !pull(tab, paste0(segment, "_AMBIG")))
        log_done()
        n_filtered <- nrow(tab)
        write_log("Rows discarded:", as.integer(n - n_filtered))
        write_log("New dimensions:", n_filtered)
    } else {
        write_log("\nNot filtering on", segment, paste0("calls", "."))
    }
}

# Write output
write_log("\nOutput dimensions:", dim(tab))
write_log("\nWriting output table to file...")
writeChangeoDb(tab, outpath)
log_done()
