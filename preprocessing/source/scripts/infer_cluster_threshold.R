#################################################################
# Compute nearest-neighbour Hamming distances for a Change-O DB #
#################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# ----------------------------------------------------------------------------
# Test for empty input by reading first line from file, then import
# ----------------------------------------------------------------------------

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
tab <- suppressMessages(readChangeoDb(inpath))
log_done()
write_log(nrow(tab), "sequence entries imported.")

# ----------------------------------------------------------------------------
# Infer clustering threshold from nearest-neighbour distribution
# ----------------------------------------------------------------------------

write_log("Inferring clustering threshold from nearest neighbours:")
models <- c("gamma-gamma", "gamma-norm", "norm-gamma", "norm-norm")
threshold_objects <- list()
plots <- list()
thresholds <- numeric(length(models))
likelihoods <- numeric(length(models))
dn <- as.numeric(tab$DIST_NEAREST)
for (n in seq(length(models))){
    model <- models[n]
    write_log("   ", paste0(model, ":"))
    threshold_objects[[n]] <- tryCatch(findThreshold(dn,
                                        method = "gmm", model = model,
                                        cutoff = cutoff),
                          error = function(e) return(e$message),
                          warning = function(w) return(w$message))
    if (!isS4(threshold_objects[[n]])) {
        thresholds[n] <- NA
        likelihoods[n] <- NA
        write_log("   ", "   ", "Thresholding failed. Error message:",
                  threshold_objects[[n]])
    } else {
        thresholds[n] <- threshold_objects[[n]]@threshold
        likelihoods[n] <- threshold_objects[[n]]@loglk
        write_log("   ", "   ", "Thresholding complete.")
        write_log("   ", "   ", "Threshold:", threshold_objects[[n]]@threshold)
        write_log("   ", "   ", "Log-likelihood:", threshold_objects[[n]]@loglk)
        write_log("   ", "   ", "p-Value:", threshold_objects[[n]]@pvalue)
    }
}
if (all(is.na(thresholds))) stop("No threshold found for any model!")

# ----------------------------------------------------------------------------
# Generate thesholding plot
# ----------------------------------------------------------------------------

write_log("\nGenerating thresholding plots:")
for (n in seq(length(models))){
    model <- models[n]
    t_o <- threshold_objects[[n]]
    t <- thresholds[n]
    l <- likelihoods[n]
    write_log("   ", paste0(model, ":"), newline = FALSE)
    if (!isS4(t_o)) {
        write_log("No threshold to plot.")
        plots[[n]] <- ggplot()
    } else {
        write_log("   ", "Generating plot...", newline=FALSE)
        # Construct plot title
        title <- paste0(model, " (thr = ", signif(t_o@threshold, 3),
                        ", loglk = ", t_o@loglk, ", p = ", 
                        signif(t_o@pvalue, 3), ")")
        plots[[n]] <- plot(t_o, binwidth = 0.005, title = title)
        if (likelihoods[n] == max(likelihoods, na.rm = TRUE)){
            # Highlight max-max-likelihood model
            plots[[n]]$layers[[1]]$aes_params$fill <- "red"
            plots[[n]]$layers[[2]]$aes_params$colour <- "purple"
            plots[[n]]$layers[[3]]$aes_params$colour <- "purple"
            plots[[n]]$layers[[4]]$aes_params$colour <- "darkgreen"
            n_max <- n
        }
        log_done()
    }
}

# Generate grid plot
write_log("Generating grid plot from solo plots...", newline=FALSE)
plot_out <- grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
                         ncol = 2, 
                         top = paste0("Distance to nearest (gmm, dist = ", 
                                      dist, ", cutoff = ", cutoff, ")"))
log_done()

# ----------------------------------------------------------------------------
# Save output files
# ----------------------------------------------------------------------------

write_log("\nWriting threshold to file...", newline = FALSE)
threshold_out <- thresholds[n_max]
write(threshold_out, outpath_threshold)
log_done()

write_log("\nWriting thresholding plot to file...", newline = FALSE)
ggsave(outpath_plot, plot = plot_out, dpi = 300, device="png")
log_done()
