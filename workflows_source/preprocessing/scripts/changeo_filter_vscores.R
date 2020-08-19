#############################################################
# Filter a Change-O DB by a V_SCORE cutoff                 ##
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
write_log(sum(as.numeric(tab$CONSCOUNT)), "surviving input reads.")

# Filter by V_SCORE
write_log("\nFiltering by V_SCORE...", newline = FALSE)
tab_filtered <- filter(tab, as.numeric(V_SCORE) >= min_vscore)
log_done()
write_log(nrow(tab_filtered),
          paste0("(", round(nrow(tab_filtered)/nrow(tab) * 100, 1),
                 "% of pre-filtered number)"),
          "sequence entries remaining.")
write_log(sum(as.numeric(tab_filtered$CONSCOUNT)),
          paste0("(", round(sum(as.numeric(tab_filtered$CONSCOUNT))/sum(as.numeric(tab$CONSCOUNT)) * 100, 1), "% of pre-filtered number)"),
          "surviving input reads.")

# Write output
write_log("\nWriting count table to file...")
write_tsv(tab_filtered, outpath)
log_done()
