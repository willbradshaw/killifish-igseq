###################################################
# Generate an entropy table from segment counts  ##
###################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")
if (!exists("unit")) unit = "log2" # Use bits by default

# Read Change-O table
write_log("\nImporting counts table...", newline = FALSE)
counts <- suppressMessages(readChangeoDb(inpath))
log_done()
write_log("Input dimensions:", dim(counts))

if (nrow(counts) == 0){
    write_log("No data; writing empty output file...")
    file.create(outpath)
    log_done()
} else {
    # Create entropy table
    write_log("\nGenerating entropy table...", newline=FALSE)
    etab <- counts %>% s_group_by(group_within, group_between) %>%
        summarise(H_DUP = entropy(DUPCOUNT, unit = "log2"),
                  H_SEQ = entropy(SEQCOUNT, unit = "log2"))
    log_done()
    write_log("Output dimensions:", dim(etab))

    # Write output
    write_log("\nWriting counts to file...")
    write_tsv(etab, outpath)
    log_done()

}

# Wildcards: group_within, group_between, set
# Params: aux
