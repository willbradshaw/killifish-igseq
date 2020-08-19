#############################################################
# Add unique sequence IDs to Change-O DB                    #
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

# Infer unique SEQUENCE_ID entries from replicate and row number
write_log("\nInferring sequence IDs...", newline = FALSE)
max_nseq <- tab %>% group_by(REPLICATE) %>% count() %>% pull(n) %>% max
max_nseq_char <- max_nseq %>% as.character %>% nchar %>% first
tab <- tab %>% group_by(REPLICATE) %>%
    mutate(SEQUENCE_ID = paste(REPLICATE,
                               str_pad(row_number(), max_nseq_char, pad="0"),
                               sep = "_")) %>% ungroup()
log_done()

# Check assigned sequence IDs are unique
write_log("\nChecking assigned IDs are unique...", newline = FALSE)
id_counts <- tab %>% group_by(SEQUENCE_ID) %>% count()
log_done()
write_log("Unique IDs: ", all(id_counts == 1))

# Write output
write_log("\nWriting count table to file...")
write_tsv(tab, outpath)
log_done()
