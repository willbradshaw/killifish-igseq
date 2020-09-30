###############################################################################
## SUPPLEMENTARY FIGURE                                                      ##
## Solo spectra of intestinal samples                                        ##
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

# Define group mapping
mapping <- list(ABX_16 = "16 weeks,\nantibiotic treated",
                SMT_16 = "16 weeks,\nsame-age transfer",
                YMT_16 = "16 weeks,\nyoung transfer",
                YI_6 = "6 weeks,\nuntreated",
                WT_16 = "16 weeks,\nuntreated",
                A_16 = "All 16-week\ngroups"
)

# Define group order and colours
groups_all <- c("YI_6", "WT_16", "ABX_16", "SMT_16", "YMT_16", "A_16")
groups_long_all <- unlist(mapping[groups_all])
palette <- c("#e78ac3", "#8da0cb", "#e5c494", "#fc8d62", "#a6d854", 
             "#F5C800")

# Configure input paths
clone_spectra_age_path   <- snakemake@input[["gut_age_clone"]]
vj_spectra_age_path      <- snakemake@input[["gut_age_vj"]]
clone_spectra_group_path <- snakemake@input[["gut_group_clone"]]
vj_spectra_group_path    <- snakemake@input[["gut_group_vj"]]

# Define output parameters
out_path <- snakemake@output[[1]]
out_width <- snakemake@params[["plot_width"]]
out_height <- snakemake@params[["plot_ratio"]] * out_width

cat("...done.\n")

#------------------------------------------------------------------------------
# IMPORT AND PROCESS DATA
#------------------------------------------------------------------------------

cat("\nPreparing data...")

# Import data
clone_spectra_age <- import_div(clone_spectra_age_path) %>% 
  mutate(CLASS = "Clonal diversity")
vj_spectra_age <- import_div(vj_spectra_age_path) %>% 
  mutate(CLASS = "VJ diversity, all clones")
clone_spectra_group <- import_div(clone_spectra_group_path) %>% 
  mutate(CLASS = "Clonal diversity")
vj_spectra_group <- import_div(vj_spectra_group_path) %>% 
  mutate(CLASS = "VJ diversity, all clones")

# Combine datasets for plotting
vj_spectra <- vj_spectra_age %>%
  mutate(GROUP = ifelse(AGE_WEEKS == 16, "A_16", "YI_6")) %>%
  bind_rows(vj_spectra_group) %>%
  mutate(GROUP_LONG = unlist(mapping[GROUP]))
clone_spectra <- clone_spectra_age %>%
  mutate(GROUP = ifelse(AGE_WEEKS == 16, "A_16", "YI_6")) %>%
  bind_rows(clone_spectra_group) %>%
  mutate(GROUP_LONG = unlist(mapping[GROUP]))
spectra <- bind_rows(clone_spectra, vj_spectra) %>%
  mutate(GROUP_LONG = factor(GROUP_LONG, groups_long_all))

cat("...done.\n")

#------------------------------------------------------------------------------
# PLOT DIVERSITY
#------------------------------------------------------------------------------

cat("\nMaking plot...")

keysize = 1.1

g_solo_spectra <-   ggplot(spectra, aes(x=Q, group=INDIVIDUAL)) +
  geom_ribbon(aes(ymin=D_LOWER_EMPIRICAL, ymax=D_UPPER_EMPIRICAL, fill=GROUP_LONG), alpha = 0.4) +
  geom_line(aes(colour=GROUP_LONG, y=D), size = 0.3) +
  scale_colour_manual(values = palette, name = "Age group (days)") +
  scale_fill_manual(values = palette, name = "Age group (days)") +
  scale_x_continuous(name="Diversity order (q)", #expand = c(0,0),
                     limits = c(0, 4),
                     breaks = seq(0, 4, 1), expand = c(0,0)) +
  scale_y_continuous(name="Diversity", expand = c(0,0),
                     limits = c(0, NA)) +
  facet_grid(CLASS~GROUP_LONG, scales = "free") +
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