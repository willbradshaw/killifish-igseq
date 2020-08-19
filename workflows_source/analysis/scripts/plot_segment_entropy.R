##############################
# Plot entropy distributions #
##############################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")
if (copy == "NULL") copy <- NULL

# Read entropy table
write_log("\nImporting entropy table...", newline = FALSE)
htab <- read_tsv(inpath)
log_done()

# Plot entropy distributions
htab[["GROUP"]] <- htab[[group_within]]
if (is.null(copy)) {
    g <- htab %>% ggplot() +
        geom_boxplot(aes(x=GROUP, y=H_SEQ, fill = factor(GROUP))) +
        theme(legend.position = "bottom") + ggtitle("SEQCOUNT")
} else if (copy == "DUPCOUNT") {
    g <- htab %>% ggplot() +
        geom_boxplot(aes(x=GROUP, y=H_DUP, fill = factor(GROUP))) +
        theme(legend.position = "bottom") + ggtitle("DUPCOUNT")
} else {
    message <- paste0("No plotting method implemented for copy value ",
                      copy, "! Choose one of (NULL,DUPCOUNT).")
    stop(message)
}
# Write output
write_log("Writing plot to image file...", newline = FALSE)
ggsave(outpath_png, plot = g, dpi = 300, device = "png")
log_done()
write_log("Saving ggplot object for editing...", newline = FALSE)
saveRDS(g, outpath_rds)
log_done()

# Inpaths: inpath
# Outpaths: outpath_png, outpath_rds
# Wildcards: set, segments, group, copy
# Params: aux
