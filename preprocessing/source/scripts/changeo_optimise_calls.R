###################################################
# Pick best V/D/J calls after germline inference ##
###################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")
if (!exists("clone_field")) clone_field <- "CLONE"
if (!exists("source_field")) source_field <- "SAMPLE"

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log(nrow(tab), "sequence entries imported.")
write_log("Input dimensions:", dim(tab))

# Optimise calls
write_log("\nOptimising V calls...", newline = FALSE)
tab <- tab %>%
    mutate(BEST_V_CALL =
           ifelse(is.na(V_CALL), GERMLINE_V_CALL,
                  ifelse(is.na(GERMLINE_V_CALL), V_CALL,
                         ifelse(str_count(GERMLINE_V_CALL, ",") <
                                str_count(V_CALL, ","),
                                GERMLINE_V_CALL, V_CALL))))
log_done()
write_log("Table dimensions:", dim(tab))
write_log("\nOptimising D calls...", newline = FALSE)
tab <- tab %>%
    mutate(BEST_D_CALL =
           ifelse(is.na(D_CALL), GERMLINE_D_CALL,
                  ifelse(is.na(GERMLINE_D_CALL), D_CALL,
                         ifelse(str_count(GERMLINE_D_CALL, ",") <
                                str_count(D_CALL, ","),
                                GERMLINE_D_CALL, D_CALL))))
log_done()
write_log("Table dimensions:", dim(tab))
write_log("\nOptimising V calls...", newline = FALSE)
tab <- tab %>%
    mutate(BEST_J_CALL =
           ifelse(is.na(J_CALL), GERMLINE_J_CALL,
                  ifelse(is.na(GERMLINE_J_CALL), J_CALL,
                         ifelse(str_count(GERMLINE_J_CALL, ",") <
                                str_count(J_CALL, ","),
                                GERMLINE_J_CALL, J_CALL))))
log_done()
write_log("Table dimensions:", dim(tab))

# Write output
write_log("\nWriting table to file...")
writeChangeoDb(tab, outpath)
log_done()
