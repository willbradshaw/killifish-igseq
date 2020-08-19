#################################################
# Generate dendrogram plots from an RDI matrix ##
#################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Read RDI matrix
write_log("\nImporting RDI matrix...", newline = FALSE)
rdi <- as.dist(read.csv(inpath, check.names = FALSE, header = TRUE,
                        row.names = 1))
log_done()
write_log("Input dimensions:", dim(rdi))

if (length(rdi) == 0){
    write_log("No data; writing empty plot files...")
    ggsave(outpath_png, plot = ggplot())
    saveRDS(ggplot(), outpath_rds)
    log_done()
} else {
    # Generate dendrogram object
    write_log("Generating dendrogram object...", newline = FALSE)
    dendrogram <- as.dendrogram(hclust(rdi))
    dendrogram_data <- dendro_data(dendrogram, type="rectangle")
    log_done()

    # Generate cluster dendrogram plot
    write_log("Generating dendrogram plot...", newline = FALSE)
    g <- ggplot() +
        geom_segment(data = segment(dendrogram_data), size = linewidth,
                     aes(x = x, y = y, xend = xend, yend = yend)) +
        scale_x_continuous(breaks = seq_along(dendrogram_data$labels$label),
                           labels = dendrogram_data$labels$label) +
        scale_y_continuous() + theme_dendro() +
        ylab("Repertoire Dissimilarity Index") +
        # TODO: Configurable xlab here?
        theme(
              axis.text.x = element_text(angle = 90, size = axis_tick_size,
                                         colour="black"),
              axis.text.y = element_text(size = axis_tick_size,
                                         colour="black"),
              axis.title.x = element_blank(),
              axis.title.y = element_text(angle = 90, size = axis_title_size,
                                          colour="black"),
              )
    log_done()

    # Write output
    write_log("Writing plot to image file...", newline = FALSE)
    ggsave(outpath_png, plot = g, dpi = 300, device = "png")
    log_done()
    write_log("Saving ggplot object for editing...", newline = FALSE)
    saveRDS(g, outpath_rds)
    log_done()
}

# Inpaths: inpath
# Outpaths: outpath_png, outpath_rds
# Wildcards: set, segments, group
# Params: linewidth, axis_tick_size, axis_title_size, aux
