###############################################################################
## SUPPLEMENTARY FIGURE                                                      ##
## Cross-replicate clone-size correlations in pilot dataset                  ##
###############################################################################

#==============================================================================
# Preamble
#==============================================================================

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Source formatting
source("scripts/aux_format-plots.R")

# Configure input paths
ctab_path    <- snakemake@input[["changeo_db"]]
# ctab_path <- "../data_processed/changeo/pilot-seqs-all.tsv.gz"

# Define output parameters
out_path <- snakemake@output[[1]]
out_width <- snakemake@params[["plot_width"]]
out_height <- snakemake@params[["plot_ratio"]] * out_width
# out_path <- "output_files/dev_si_pilot-cor.png"
# out_width <- 11
# out_ratio <- 1.5
# out_height <- out_width * out_ratio

#==============================================================================
# Auxiliary functions
#==============================================================================

cat("done.\n")

#==============================================================================
# Import data
#==============================================================================

cat("\nImporting data...")
ctab <- suppressMessages(read_tsv(ctab_path))
cat("done.\n")

#==============================================================================
# Process data
#==============================================================================

cat("\nProcessing data...")

# Count clones per individual
tab_cl <- ctab %>% mutate(REP = sub("\\d-\\d\\d", "", REPLICATE)) %>%
  filter(!is.na(CLONE)) %>%
  group_by(INDIVIDUAL, REP, CLONE) %>% 
  summarise(CLNCOUNT = n(), DUPCOUNT = sum(DUPCOUNT), CONSCOUNT = sum(CONSCOUNT))

tab_cl_spread <- spread(tab_cl, REP, CLNCOUNT) %>%
  mutate(bio = ifelse(is.na(bio), 0, bio),
         lib = ifelse(is.na(lib), 0, lib),
         orig = ifelse(is.na(orig), 0, orig)) %>%
  group_by(INDIVIDUAL, CLONE) %>%
  summarise(DUPCOUNT = sum(DUPCOUNT), CONSCOUNT = sum(CONSCOUNT),
            bio = sum(bio), lib = sum(lib), orig = sum(orig)) %>%
  mutate(N_ABSENT = ((bio==0) + (lib==0) + (orig == 0)),
         N_PRESENT = 3-N_ABSENT,
         CLNCOUNT = (bio+lib+orig),
         CLNCOUNT_AVG = (bio+lib+orig)/(N_PRESENT))

# Get inter-replicate correlation in clone sizes
tab_cl_cor <- tab_cl_spread %>%
  summarise(BL = cor(bio, lib), BO = cor(bio, orig), LO = cor(lib, orig)) %>%
  melt(id.vars = "INDIVIDUAL", variable.name = "COMPARISON", 
       value.name = "R")

cat("done.\n")

#==============================================================================
# Make plots 
#==============================================================================

cat("\nGenerating boxplots...")

# Plot distribution of inter-replicate correlation coefficients
g_cor <- ggplot(tab_cl_cor, aes(x=COMPARISON, y=R)) + 
  geom_boxplot(aes(fill = COMPARISON), outlier.shape = NA) +
  geom_point(size = 1, alpha = 0.4, shape = 16) +
  scale_x_discrete(labels = c("1 vs 2", "1 vs 3", "2 vs 3"),
                   name = "Inter-replicate comparison") + 
  scale_y_continuous(breaks = seq(0,1,0.1),
                     limits = c(0.6,1), name = "Correlation coefficient (r)") + 
  scale_fill_brewer(palette = "Set1") +
  theme_base +
  theme(legend.position = "none",
        panel.grid = element_line(colour = "#F5F5F5"),
        aspect.ratio = 0.5,
        )

cat("done.\n")

cat("\nGenerating scatter plots...")

g_inter <- ggplot(tab_cl_spread) + 
  geom_point(aes(x=orig, y=bio, colour = INDIVIDUAL),
             size = 1, alpha = 0.7, shape = 16) +
  facet_wrap(~INDIVIDUAL, scales = "free") +
  scale_colour_manual(values = pilot_palette, name = "Individual") +
  scale_x_log10(name = "Replicate 1", limits = c(1,1000)) +
  scale_y_log10(name = "Replicate 2", limits = c(1,1000)) +
  theme_base + theme(
    legend.position = "none",
    aspect.ratio = 1,
    strip.text.x = element_text(face="bold"),
    )

cat("done.\n")

cat("\nAssembling output grid...")
fig_out <- plot_grid(g_cor, g_inter, ncol = 1, nrow = 2, labels = "AUTO",
                        label_size = fontsize_base * fontscale_label,
                        label_fontfamily = titlefont, label_colour = "black",
                        label_fontface = "plain", vjust=1.1,
                        rel_heights = c(0.5,1))
cat("done.\n")

#==============================================================================
# SAVE PLOT
#==============================================================================

cat("\nSaving output...")
save_fig(out_path, fig_out, out_height, out_width)
cat("...done.\n")