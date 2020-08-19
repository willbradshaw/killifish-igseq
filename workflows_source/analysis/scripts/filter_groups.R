###################################################################
# Filter a Change-O DB by group specifications                    #
###################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Specify default values
if (!exists("keep_complement")) keep_complement <- TRUE

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Filter columns
write_log("\nFiltering on specified columns...")
for (col in names(filter_groups)){
    write_log("   ", col)
    vals_all <- unique(tab[[col]])
    vals_listed <- filter_groups[[col]]
    keep <- vals_all %in% vals_listed
    if (keep_complement) keep <- !keep
    vals_keep <- vals_all[keep]
    tab <- tab[tab[[col]] %in% vals_keep,]
}
log_done()

# Write output
write_log("\nOutput dimensions:", dim(tab))
write_log("\nWriting output table to file...", newline = FALSE)
writeChangeoDb(tab, outpath)
log_done()
