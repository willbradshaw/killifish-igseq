#################################################################
# Find best thresholding model and make plots                  ##
#################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# ----------------------------------------------------------------------------
# Import threshold objects
# ----------------------------------------------------------------------------

models <- c("gamma-gamma", "gamma-norm", "norm-gamma", "norm-norm")
write_log("\nImporting threshold objects...", newline = FALSE)
threshold_objects <- lapply(inpaths, readRDS)
names(threshold_objects) <- models
log_done()

write_log("Discarding failed models...", newline = FALSE)
threshold_objects <- threshold_objects[sapply(threshold_objects, isS4)]
log_done()
write_log(length(threshold_objects), "models remaining.")

# ----------------------------------------------------------------------------
# Extract thresholds and log-likelihoods
# ----------------------------------------------------------------------------

thresholds <- sapply(threshold_objects, function(t) t@threshold)
likelihoods <- sapply(threshold_objects, function(t) t@loglk)

# ----------------------------------------------------------------------------
# Generate thesholding plots
# ----------------------------------------------------------------------------

write_log("\nGenerating thresholding plots:")
plots <- list()
for (n in seq(length(models))){
    model <- models[n]
    write_log("   ", paste0(model, ":"), newline = FALSE)
    if (!model %in% names(threshold_objects)){
        write_log("No threshold to plot.")
        plots[[n]] <- ggplot()
    } else {
        t_o <- threshold_objects[[n]]
        t <- thresholds[n]
        l <- likelihoods[n]
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
ggsave(outpath_plot, plot = plot_out, dpi = 320, device="png")
log_done()
