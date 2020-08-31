###############################################################################
## SUPPLEMENTARY FIGURE                                                      ##
## KW p-values of whole-body ageing samples                                  ##
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
clone_kw_path    <- snakemake@input[["clone"]]
vj_all_kw_path   <- snakemake@input[["vj_all"]]
vj_large_kw_path <- snakemake@input[["vj_large"]]
vj_small_kw_path <- snakemake@input[["vj_small"]]

# Define output parameters
out_path <- snakemake@output[[1]]
out_width <- snakemake@params[["plot_width"]]
out_height <- snakemake@params[["plot_ratio"]] * out_width

# Modify theme
theme_base <- theme_base + theme(
  legend.text = element_text(margin = margin(r = 0.2, unit = "cm"))
)

cat("...done.\n")

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

cat("\nPreparing data...")

clone_kw <- import_div(clone_kw_path) %>% 
  mutate(CLASS = "Clonal diversity")
vj_kw_all <- import_div(vj_all_kw_path) %>% 
  mutate(CLASS = "VJ diversity,\nall clones")
vj_kw_small <- import_div(vj_small_kw_path) %>% 
  mutate(CLASS = "VJ diversity,\nsmall clones")
vj_kw_large <- import_div(vj_large_kw_path) %>% 
  mutate(CLASS = "VJ diversity,\nlarge clones")

kw <- bind_rows(clone_kw, vj_kw_all, vj_kw_small,
                     vj_kw_large) %>%
  filter(METHOD == "KW permutation test")

cat("...done.\n")

#------------------------------------------------------------------------------
# PLOT P-VALUES
#------------------------------------------------------------------------------

cat("\nPreparing plot...")

keysize = 0.9

g_kw <- ggplot(kw, aes(x=Q, y=P, colour = CLASS)) +
  geom_hline(yintercept = c(0.05, 0.01, 0.001), linetype = "dotted", 
             colour = "red", size = 0.5) +
  geom_line() +
  scale_y_log10(name = "P-value (Kruskal-Wallis permutation test)",
                breaks = c(0.05, 0.01, 0.001, 0.3, 1)) +
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