###############################################################################
## SUPPLEMENTARY FIGURE                                                      ##
## KW p-values of intestinal samples                                         ##
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
clone_kw_age_path <- snakemake@input[["gut_age_clone"]]
vj_kw_age_path <- snakemake@input[["gut_age_vj"]]
clone_kw_group_path <- snakemake@input[["gut_group_clone"]]
vj_kw_group_path <- snakemake@input[["gut_group_vj"]]

# Define output parameters
out_path <- snakemake@output[[1]]
out_width <- snakemake@params[["plot_width"]]
out_height <- snakemake@params[["plot_ratio"]] * out_width

cat("...done.\n")

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

cat("\nPreparing data...")

clone_kw_age <- import_div(clone_kw_age_path) %>% 
  mutate(CLASS = "Clonal diversity", TEST = "Age effect")
vj_kw_age <- import_div(vj_kw_age_path) %>% 
  mutate(CLASS = "VJ diversity, all clones", TEST = "Age effect")
clone_kw_group <- import_div(clone_kw_group_path) %>% 
  mutate(CLASS = "Clonal diversity", TEST = "Treatment effect")
vj_kw_group <- import_div(vj_kw_group_path) %>% 
  mutate(CLASS = "VJ diversity, all clones", TEST = "Treatment effect")

kw <- bind_rows(clone_kw_age, vj_kw_age, clone_kw_group, vj_kw_group) %>%
  filter(METHOD == "KW permutation test") %>%
  mutate(TEST = factor(TEST, levels = c("Treatment effect", "Age effect")))

cat("...done.\n")

#------------------------------------------------------------------------------
# PLOT P-VALUES
#------------------------------------------------------------------------------

cat("\nMaking plot...")

keysize = 0.9

g_kw <- ggplot(kw, aes(x=Q, y=P, colour = CLASS)) +
  geom_hline(yintercept = c(0.05, 0.01, 0.001), linetype = "dotted", 
             colour = "red", size = 0.5) +
  geom_line() +
  scale_y_log10(name = "P-value (Kruskal-Wallis permutation test)",
                breaks = c(0.05, 0.01, 0.001, 0.3, 1)) +
  facet_grid(.~TEST) +
  scale_colour_brewer(palette = "Set1", name = "Diversity type") +
  coord_fixed() +
  guides(colour = guide_legend(keyheight=keysize, keywidth = keysize,
                               nrow = 1)) +
  theme_base

cat("...done.\n")

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

cat("\nSaving output...")

save_fig(out_path, g_kw, out_height, out_width)

cat("...done.\n")
