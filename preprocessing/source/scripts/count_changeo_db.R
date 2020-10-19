#############################################################
# Count sequence number and read usage for a Change-O table #
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

# If no CONSCOUNT column, count rows instead
if (! "CONSCOUNT" %in% colnames(tab)) tab$CONSCOUNT <- 1
if (! "NREADS_RAW" %in% colnames(tab)) tab$NREADS_RAW <- NA

# Fix READS_RAW column (for multiple values)
tab <- mutate(tab, NREADS_RAW = str_split(NREADS_RAW, ","))

# Make count table
write_log("\nGenerating count table...", newline = FALSE)
count_tab <- tab %>% s_group_by(group_field) %>%
    summarise(NSEQS = n(), NREADS = sum(as.integer(CONSCOUNT)),
              NREADS_RAW = sum(unique(as.integer(unlist(NREADS_RAW))))) %>%
    mutate(STAGE = stage, CLUSTER_BARCODES = cluster_barcodes,
           CLUSTER_SETS = cluster_sets,
           NREADS_PC = NREADS/NREADS_RAW) # Fields from params
log_done()

# Write output
write_log("Writing count table to file...")
write_tsv(count_tab, outpath)
log_done()
