######################################################
# Generate a segment count table from a Change-O DB ##
######################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
tab <- suppressMessages(readChangeoDb(inpath)) %>% mutate(SEQCOUNT = 1)
log_done()
write_log("Input dimensions:", dim(tab))

if (nrow(tab) == 0){
    write_log("No data; writing empty output file...")
    file.create(outpath)
    log_done()
} else {

    # Determine formulae
    f1_suffix <- paste0(group_within, "+",group_between, "+", call_field)
    f1_dup <- as.formula(paste0("DUPCOUNT~", f1_suffix))
    f1_seq <- as.formula(paste0("SEQCOUNT~", f1_suffix))
    f2 <- as.formula(paste0(group_within, "+", group_between, "~", call_field))

    # Compute frequencies from input
    write_log("\nComputing count frequencies...", newline=FALSE)
    counts_dup <- tab %>% aggregate(f1_dup, ., sum) %>%
        dcast(f2, fun.aggregate = sum, value.var = "DUPCOUNT") %>%
        melt(id.vars = c(group_within, group_between),
             variable.name = call_field, value.name = "DUPCOUNT")
    counts_seq <- tab %>% aggregate(f1_seq, ., sum) %>%
        dcast(f2, fun.aggregate = sum, value.var = "SEQCOUNT") %>%
        melt(id.vars = c(group_within, group_between),
             variable.name = call_field, value.name = "SEQCOUNT")
    counts <- inner_join(counts_dup, counts_seq, 
                         by = c(group_within, group_between, call_field))
    counts[[call_field]] <- as.character(counts[[call_field]])
    log_done()
    write_log("Output dimensions:", dim(counts))

    # Write output
    write_log("\nWriting counts to file...")
    write_tsv(counts, outpath)
    log_done()

}

# Wildcards: group_within, group_between, set
# Params: aux, call_field
