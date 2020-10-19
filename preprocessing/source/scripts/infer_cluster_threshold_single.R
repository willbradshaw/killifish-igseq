#################################################################
# Compute clonotype distance threshold under a single model     #
#################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")
write_log("Distance metric:", dist)
write_log("Fitting model:", model)
write_log("Cutoff method:", cutoff)

# ----------------------------------------------------------------------------
# Import Change-O DB
# ----------------------------------------------------------------------------

write_log("\nImporting Change-O table...", newline = FALSE)
tab <- suppressMessages(readChangeoDb(inpath))
log_done()
write_log(nrow(tab), "sequence entries imported.")

# ----------------------------------------------------------------------------
# Infer clustering threshold from nearest-neighbour distribution
# ----------------------------------------------------------------------------
dn <- as.numeric(tab$DIST_NEAREST)

write_log("Inferring clustering threshold from nearest neighbours:")
write_log("   ", paste0(model, ":"))
threshold_object <- tryCatch(findThreshold(dn,
                                    method = "gmm", model = model,
                                    cutoff = cutoff),
                      error = function(e) return(e$message),
                      warning = function(w) return(w$message))
if (!isS4(threshold_object)) {
    write_log("   ", "   ", "Thresholding failed. Error message:",
              threshold_object)
} else {
    threshold <- threshold_object@threshold
    likelihood <- threshold_object@loglk
    write_log("   ", "   ", "Thresholding complete.")
    write_log("   ", "   ", "Threshold:", threshold_object@threshold)
    write_log("   ", "   ", "Log-likelihood:", threshold_object@loglk)
    write_log("   ", "   ", "p-Value:", threshold_object@pvalue)
}

# ----------------------------------------------------------------------------
# Serialise and save threshold object
# ----------------------------------------------------------------------------

write_log("\nWriting threshold object to file...", newline = FALSE)
saveRDS(threshold_object, outpath)
log_done()
