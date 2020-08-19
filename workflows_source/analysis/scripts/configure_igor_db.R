######################################################################
# Prepare sequence DB for sequence extraction for IGoR analysis      #
######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")
if (!exists("collapsed")) collapsed <- FALSE

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(MU_COUNT_CDR_R = "d", MU_COUNT_CDR_S = "d",
            MU_COUNT_FWR_R = "d", MU_COUNT_FWR_S = "d",
            MU_EXPECTED_CDR_R = "d", MU_EXPECTED_CDR_S = "d",
            MU_EXPECTED_FWR_R = "d", MU_EXPECTED_FWR_S = "d",
            .default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Filter individuals
write_log("\nFiltering individuals...", newline = FALSE)
tab <- tab %>% filter(INDIVIDUAL %in% as.character(individuals))
log_done()
write_log(nrow(tab), "sequence entries remaining.")

# Filter groups
write_log("\nFiltering", group_field, "groups...", newline = FALSE)
tab <- tab %>% filter(!!as.name(group_field) %in% as.character(groups))
log_done()
write_log(nrow(tab), "sequence entries remaining.")

# Remove gap characters from clonal sequence
write_log("\nStripping ambiguous/gap characters from sequences...",
          newline = FALSE)
if (collapsed){
    tab <- tab %>% mutate(SEQUENCE_OUT = gsub("[^ACGTN]", "", CLONAL_SEQUENCE))
} else {
    tab <- tab %>% mutate(SEQUENCE_OUT = gsub("[^ACGTN]", "", SEQUENCE_IMGT))
}
log_done()
write_log(nrow(tab), "sequence entries remaining.")

# Write output
write_log("\nOutput dimensions:", dim(tab))
write_log("\nWriting output table to file...", newline = FALSE)
writeChangeoDb(tab, outpath)
log_done()
