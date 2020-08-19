###################################################################
# Filter a Change-O DB by clone size                              #
###################################################################
# NB: Currently only groups by individual; needs refinement if
# want to group by replicate instead / as well

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Specify default values
if (!exists("min_size")) min_size <- 1
if (!exists("max_size")) max_size <- Inf
if (!exists("keep_complement")) keep_complement <- TRUE

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Count unique sequences per clone per individual
write_log("\nCounting clone sizes...")
clone_sizes <- tab %>% group_by(INDIVIDUAL, CLONE, SEQUENCE_IMGT) %>%
    summarise(CLNCOUNT = 1) %>% group_by(INDIVIDUAL, CLONE) %>%
    summarise(CLNCOUNT = sum(CLNCOUNT))
tab_cs <- full_join(tab, clone_sizes, by = c("INDIVIDUAL", "CLONE"))
log_done()
write_log("Dimensions:", dim(tab_cs))

# Filter sequences with missing clonal identities
write_log("\nDiscarding sequences with missing clonal identities...")
tab_nona <- tab_cs %>% filter(!is.na(CLONE))
log_done()
write_log("Dimensions:", dim(tab_nona))
write_log("Entries discarded:", nrow(tab_cs) - nrow(tab_nona))

# Filter sequences by clone size
write_log("\nDiscarding sequences by clone size...")
tab_sized <- tab_nona %>% filter(between(CLNCOUNT, min_size, max_size))
log_done()
write_log("Dimensions:", dim(tab_sized))
write_log("Entries discarded:", nrow(tab_nona) - nrow(tab_sized))

# Write output
write_log("\nWriting output table to file...", newline = FALSE)
writeChangeoDb(tab_sized, outpath)
log_done()
