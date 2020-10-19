################################################################
# Add VJ and VDJ combined segment call fields to a Change-O DB #
################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(HAS_V = "l", HAS_D = "l", HAS_J = "l",
            V_AMBIG = "l", D_AMBIG = "l", J_AMBIG = "l",
            .default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Add VJ combined calls
write_log("\nAdding combined VJ segment calls...", newline = FALSE)
tab_vj <- tab %>%
    unite(VJ_CALL, V_CALL, J_CALL, sep="/", remove=FALSE) %>%
    unite(GERMLINE_VJ_CALL, GERMLINE_V_CALL, GERMLINE_J_CALL,
          sep="/", remove=FALSE) %>%
    unite(BEST_VJ_CALL, BEST_V_CALL, BEST_J_CALL,
          sep="/", remove=FALSE) %>%
    mutate(HAS_VJ = HAS_V & HAS_J, VJ_AMBIG = V_AMBIG | J_AMBIG)
log_done()
write_log("New dimensions:", dim(tab_vj))

# Add VDJ combined calls
write_log("\nAdding combined VDJ segment calls...", newline = FALSE)
tab_vdj <- tab_vj %>%
    unite(VDJ_CALL, V_CALL, D_CALL, J_CALL, sep="/", remove=FALSE) %>%
    unite(GERMLINE_VDJ_CALL, GERMLINE_V_CALL, GERMLINE_D_CALL,
          GERMLINE_J_CALL, sep="/", remove=FALSE) %>%
    unite(BEST_VDJ_CALL, BEST_V_CALL, BEST_D_CALL, BEST_J_CALL,
          sep="/", remove=FALSE) %>%
    mutate(HAS_VDJ = HAS_VJ & HAS_D, VDJ_AMBIG = VJ_AMBIG | D_AMBIG)
log_done()
write_log("New dimensions:", dim(tab_vdj))

# Write output
write_log("\nWriting output table to file...")
writeChangeoDb(tab_vdj, outpath)
log_done()
