#############################################################
# Count sequence number and read usage for a Change-O table #
#############################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")
if (!exists("collapsed")) collapsed <- TRUE
write_log("Counting reads from", 
          ifelse(collapsed, "provided CONSCOUNT values.",
                 "number of input rows."))
if (!exists("count_raw")) count_raw <- !exists("inpath_raw")
write_log("Computing raw read count from",
          ifelse(count_raw, "number of input rows.",
                 "provided raw-read count table."))

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
tab <- suppressMessages(readChangeoDb(inpath_tab))
log_done()
write_log(nrow(tab), "sequence entries imported.")

# Read raw read count from count table, if necessary
write_log("\nObtaining raw read count...", newline = FALSE)
if (count_raw){
    raw <- nrow(tab)
} else {
    raw <- read_tsv(inpath_raw, col_names = FALSE)[1,2] %>%
        as.numeric()
}
log_done()
write_log("Raw read count:", raw)

# Make count table
write_log("\nGenerating count table...", newline = FALSE)
count_tab <- tibble(
                    STAGE = stage,
                    SEQCOUNT = nrow(tab),
                    CONSCOUNT = ifelse(collapsed, sum(tab$CONSCOUNT), nrow(tab)),
                    SAMPLE = sample, # From params
                    SIZE = size, # From params
                    ITERATION = iter, # From params
                    CLUSTER_BARCODES = cluster_barcodes, # From params
                    CLUSTER_SETS = cluster_sets, # From params
                    CONSCOUNT_RAW = raw,
                    ) %>% mutate(CONSCOUNT_PC = CONSCOUNT/CONSCOUNT_RAW)
log_done()

# Write output
write_log("Writing count table to file...")
write_tsv(count_tab, outpath)
log_done()
