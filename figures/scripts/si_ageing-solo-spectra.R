###############################################################################
## SUPPLEMENTARY FIGURE                                                      ##
## Solo spectra of whole-body ageing samples                                 ##
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
clone_spectra_path   <- snakemake@input[["clone"]]
vj_all_spectra_path      <- snakemake@input[["vj_all"]]
vj_large_spectra_path      <- snakemake@input[["vj_large"]]
vj_small_spectra_path      <- snakemake@input[["vj_small"]]

# Define output parameters
out_path <- snakemake@output[[1]]
out_width <- snakemake@params[["plot_width"]]
out_height <- snakemake@params[["plot_ratio"]] * out_width

cat("...done.\n")

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

cat("\nPreparing data...")

clone_spectra <- import_div(clone_spectra_path) %>% 
  mutate(CLASS = "Clonal diversity")
vj_spectra_all <- import_div(vj_all_spectra_path) %>% 
  mutate(CLASS = "VJ diversity,\nall clones")
vj_spectra_small <- import_div(vj_small_spectra_path) %>% 
  mutate(CLASS = "VJ diversity,\nsmall clones")
vj_spectra_large <- import_div(vj_large_spectra_path) %>% 
  mutate(CLASS = "VJ diversity,\nlarge clones")

spectra <- bind_rows(clone_spectra, vj_spectra_all, vj_spectra_small,
                     vj_spectra_large) %>%
  mutate(AGE_DAYS = factor(AGE_DAYS, levels = ageing_groups))

cat("...done.\n")

#------------------------------------------------------------------------------
# PLOT DIVERSITY
#------------------------------------------------------------------------------

cat("\nGenerating plot...")

keysize = 0.7

g_solo_spectra <-   ggplot(spectra, aes(x=Q, group=INDIVIDUAL)) +
  geom_ribbon(aes(ymin=D_LOWER_EMPIRICAL, 
                  ymax=D_UPPER_EMPIRICAL, fill=AGE_DAYS), alpha = 0.4) +
  geom_line(aes(colour=AGE_DAYS, y=D), size = 0.3) +
  scale_colour_manual(values = ageing_palette, name = "Age group (days)") +
  scale_fill_manual(values = ageing_palette, name = "Age group (days)") +
  scale_x_continuous(name="Diversity order (q)", #expand = c(0,0),
                     limits = c(0, 4),
                     breaks = seq(0, 4, 1), expand = c(0,0)) +
  scale_y_continuous(name="Diversity", expand = c(0,0),
                     limits = c(0, NA)) +
  facet_grid(CLASS~AGE_DAYS, scales = "free",
             labeller = labeller(AGE_DAYS = function(x) paste(x, "days"))) +
  guides(colour = guide_legend(keyheight=keysize, keywidth = keysize),
         fill   = guide_legend(keyheight=keysize, keywidth = keysize)) +
  theme_base + theme(
    legend.position = "bottom",
    panel.grid = element_line(colour = "grey92")
    )

cat("...done.\n")

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

cat("\nSaving output...")
save_fig(out_path, g_solo_spectra, out_height, out_width)
cat("...done.\n")