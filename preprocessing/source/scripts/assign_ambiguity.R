######################################################################
# Label rows of a Change-O DB with ambiguity status of segment calls #
######################################################################

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

# Test assigned V/D/J segments
write_log("\nCounting assigned segments and testing for ambiguity...", 
          newline = FALSE)
tab <- mutate(tab,
              HAS_V = !is.na(BEST_V_CALL), HAS_D = !is.na(BEST_D_CALL),
              HAS_J = !is.na(BEST_J_CALL),
              N_V = str_count(BEST_V_CALL, ",") + 1,
              N_D = str_count(BEST_D_CALL, ",") + 1,
              N_J = str_count(BEST_J_CALL, ",") + 1,
              V_AMBIG = N_V > 1, D_AMBIG = N_D > 1, J_AMBIG = N_J > 1,
              HAS_CLONE = !is.na(CLONE))
log_done()
write_log("\nProportion of sequences with valid V call:", mean(tab$HAS_V))
write_log("Proportion of sequences with valid D call:", mean(tab$HAS_D))
write_log("Proportion of sequences with valid J call:", mean(tab$HAS_J))
write_log("Proportion of sequences with unambiguous V call:",
          mean(tab$V_AMBIG, na.rm = TRUE))
write_log("Proportion of sequences with unambiguous D call:",
          mean(tab$D_AMBIG, na.rm = TRUE))
write_log("Proportion of sequences with unambiguous J call:",
          mean(tab$J_AMBIG, na.rm = TRUE))
write_log("Average number of (valid) V candidates:", mean(tab$N_V, na.rm = TRUE))
write_log("Average number of (valid) D candidates:", mean(tab$N_D, na.rm = TRUE))
write_log("Average number of (valid) J candidates:", mean(tab$N_J, na.rm = TRUE))

# Write output
write_log("\nExporting", nrow(tab), "sequence entries...")
write_tsv(tab, outpath)
log_done()
