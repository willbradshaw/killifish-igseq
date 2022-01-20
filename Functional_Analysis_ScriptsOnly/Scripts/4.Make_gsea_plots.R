#==============================================================================
# Preamble
#==============================================================================

# Import libraries
library(tidyverse)
library(patchwork)

# Specify parameters
q_values_raw <- c(0, 1, 1.5, 2, 3, 4)
q_values_str <- format(q_values_raw, nsmall=2)
immune_cols <- c("immune_system_processes", "regulation_of_immune_system_processes",
                 "immune_response_regulating_signaling_pathway",
                 "cytokine_production",
                 "cytokine_mediated_signaling_pathway",
                 "defense_response", "regulation_of_defense_response",
                 "leukocyte_proliferation")

n_top <- 10 # Number of top GO terms to display for each diversity order
plot_width <- 16
plot_ratio_top <- 2/3
plot_ratio_all <- 2/3
top_buffer_x <- 0.05

# Specify input paths
in_paths <- sapply(q_values_str, function(q)
  paste0("Output/GSEA_analysis/gsea_simplified_ImmuneTerms_Q", q, ".csv"))

# Specify output paths
out_path_robust <- "Figures/gut_gsea_robust.csv"
out_path_top_pos <- "Figures/gut_gsea_top_pos.png"
out_path_top_neg <- "Figures/gut_gsea_top_neg.png"
out_path_all_pos <- "Figures/gut_gsea_all_pos.png"
out_path_all_neg <- "Figures/gut_gsea_all_neg.png"
#------------------------------------------------------------------------------
# Fonts
#------------------------------------------------------------------------------

font <- "sans" # Main figure font
titlefont <- font # Font for axis titles etc. (if different)
fontsize_base <- 6 # Basic figure font size
fontscale_title <- 1 # Default axis-title scale relative to regular font
fontscale_main <- 1 # Default plot-title scale
fontscale_label <- 2 # Default subfigure-label scale (A, B, etc)
fontscale_legend <- 1

#------------------------------------------------------------------------------
# Themes
#------------------------------------------------------------------------------

axis_width <- 0.3
theme_base <-   theme_bw() + theme(
  legend.position = "bottom",
  axis.text = element_text(size = fontsize_base, family = font, 
                           colour = "black"),
  axis.title.y = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(r=1,unit="mm")),
  axis.title.x = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(t=1,unit="mm")),
  legend.text = element_text(size = fontsize_base * fontscale_legend, 
                             family = font, colour = "black"),
  legend.title = element_text(size = fontsize_base * fontscale_legend,
                              family = font, colour = "black",
                              face = "bold", vjust=0.5, 
                              margin=margin(r=3, unit="mm")),
  panel.grid = element_blank(),
  plot.margin = margin(t=0.5, l=0.2, r=0.4, b = 0.2, unit="cm"),
  strip.background = element_blank(),
  legend.justification = "center",
  panel.border = element_blank(),
  axis.line.x = element_line(size=axis_width, colour="black"),
  axis.line.y = element_line(size=axis_width, colour="black"),
  axis.ticks.x = element_line(size=axis_width, colour="black"),
  axis.ticks.y = element_line(size=axis_width, colour="black"),
  strip.text.x = element_text(size = fontsize_base, family = font, colour = "black"),
  strip.text.y = element_text(size = fontsize_base, family = font, colour = "black"),
)

#==============================================================================
# Import and collate data
#==============================================================================

# Import data
n_sets <- length(in_paths)
gsea_data_raw <- lapply(in_paths, function(x) read_csv(x, col_types = cols()))

# Annotate with diversity orders
n_sets <- length(in_paths)
q_values_lab <- format(q_values_raw, nsmall=1) %>% fct_inorder()
gsea_data_div <- lapply(1:n_sets, function(n) 
  gsea_data_raw[[n]] %>% mutate(diversity_order = q_values_lab[n]))

# Collate & filter columns
fix_columns <- function(tab){
  tab %>% select(-1, -X, -setSize, -enrichmentScore, -rank, -leading_edge,
                 -core_enrichment) %>% # Drop unused columns
    mutate(go_number = as.numeric(sub("GO:", "", ID))) # Extract GO term numbers (for arbitrary tiebreaking))
}
gsea_data <- bind_rows(gsea_data_div) %>% fix_columns

#==============================================================================
# Identify all immune terms
#==============================================================================

collate_immune <- function(tab, cols = immune_cols){
  is_immune_vec <- rep(FALSE, nrow(tab))
  for (col in cols){
    is_immune_vec <- is_immune_vec | tab[[col]]
  }
  tab_immune <- mutate(tab, is_immune = is_immune_vec)
  for (col in cols){tab_immune[[col]] <- NULL}
  return(tab_immune)
}

gsea_data_immune <- collate_immune(gsea_data)

#==============================================================================
# Separate positive and negative terms and rank by NES
#==============================================================================

rank_nes <- function(tab){
  # Rank by abs(NES) regardless of sign
  mutate(tab, NES_abs = abs(NES), NES_pos = NES >= 0) %>% 
    arrange(desc(NES_abs), go_number) %>% # Sort by absolute NES in descending order
    mutate(NES_rank_all = row_number()) %>%
    group_by(NES_pos, .add = TRUE) %>%
    mutate(NES_rank_sign = row_number()) %>%
    ungroup
}
gsea_data_nes <- gsea_data_immune %>% group_by(diversity_order) %>% rank_nes

#==============================================================================
# Identify robustly enriched terms (in diversity-order data)
#==============================================================================

# Group by GO ID and count diversity orders
gsea_data_ndiv <- gsea_data_nes %>%
  group_by(ID, Description, is_immune) %>% 
  summarise(n_diversity_orders = n(), 
            avg_rank = mean(NES_rank_all),
            .groups = "drop_last")

# Filter by immune status and number of orders
gsea_data_fdiv <- gsea_data_ndiv %>%
  filter(is_immune, n_diversity_orders >= 4) %>%
  arrange(desc(n_diversity_orders), avg_rank)

# Extract information about all robustly enriched immune terms
gsea_data_robust_full <- lapply(gsea_data_fdiv$ID, function(id)
  gsea_data_nes %>% ungroup %>% filter(ID == id) %>% 
    select(ID, Description, diversity_order, NES, NES_rank_all, 
           NES_rank_sign, p.adjust)
)

# Summarise for display
msd <- function(x, n=0){
  paste(round(mean(x), n), "±", round(sd(x), n+1))
}
gsea_data_robust_summ <- gsea_data_robust_full %>% lapply(
  function(tab) tab %>% summarise(
    ID = ID[1], description = Description[1], 
    diversity_orders = paste(diversity_order, collapse = ","),
    n_orders = n(),
    NES = msd(NES, 1), rank_sign = msd(NES_rank_sign),
    rank_all = msd(NES_rank_all), p_adjust = msd(p.adjust, 2))
) %>% bind_rows

#==============================================================================
# Plot all significant terms (diversity-order data)
#==============================================================================

# Theme
theme_all <- theme_base + theme(
  axis.title.y = element_text(margin = margin(r=1.5, unit = "mm")),
  panel.grid.major = element_line(colour = "grey92", size = rel(0.5)), # Inherits from line
  legend.title = element_blank(),
  legend.box.margin = margin(0),
  legend.margin = margin(0),
  legend.box.spacing = unit(1, "mm")
)

# Prepare data
gsea_data_all <- gsea_data_nes %>%
  mutate(immune_lab = ifelse(is_immune, "Immune terms", "Non-immune terms"),
         immune_lab = factor(immune_lab, levels = rev(c("Immune terms", "Non-immune terms"))))
gsea_data_all_pos <- filter(gsea_data_all, NES_pos)
gsea_data_all_neg <- filter(gsea_data_all, !NES_pos)

# Define palette
palette_immune <- c("#BBBBBB", "#1f78b4")
scale_colour_immune <- partial(scale_colour_manual, name = NULL,
                               values = palette_immune)

# Positively enriched
g_all_pos <- ggplot(gsea_data_all_pos, aes(x=NES_rank_sign, y=NES)) +
  geom_vline(xintercept = 10, linetype = "dotted", colour = "red", size = 0.5) +
  geom_line(aes(group = diversity_order), colour = last(palette_immune)) +
  geom_point(aes(colour = immune_lab)) +
  facet_grid(diversity_order ~ ., 
             labeller = labeller(diversity_order = as_labeller(function(x) paste0("q = ", x)))) +
  scale_x_continuous(name = "Rank", breaks = seq(0,1000,10), limits = c(0,100)) + 
  scale_y_continuous(name = "Normalised Enrichment Score") +
  scale_colour_immune() +
  theme_all

# Negatively enriched
g_all_neg <- ggplot(gsea_data_all_neg, aes(x=NES_rank_sign, y=NES)) +
  geom_vline(xintercept = 10, linetype = "dotted", colour = "red", size = 0.5) +
  geom_line(aes(group = diversity_order), colour = last(palette_immune)) +
  geom_point(aes(colour = immune_lab)) +
  facet_grid(diversity_order ~ ., 
             labeller = labeller(diversity_order = as_labeller(function(x) paste0("q = ", x)))) +
  scale_x_continuous(name = "Rank", breaks = seq(0,1000,10), limits = c(0,100)) + 
  scale_y_continuous(name = "Normalised Enrichment Score") +
  scale_colour_immune() +
  theme_all

#==============================================================================
# Top-N plots (diversity-order data)
#==============================================================================

# Prepare theme
theme_stack <- theme_base + theme(
  axis.title.y = element_blank(),
  panel.grid.major = element_line(colour = "grey92", size = rel(0.5)), # Inherits from line
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                             margin = margin(t=0.3, b = 0.3, unit="mm")),
  legend.title = element_blank(),
  legend.box.margin = margin(0),
  legend.margin = margin(0),
  legend.box.spacing = unit(1, "mm")
)

# Prepare data
gsea_data_top <- filter(gsea_data_all, NES_rank_sign <= n_top)
gsea_data_top_pos <- filter(gsea_data_top, NES_pos) %>% ungroup %>%
  mutate(label = fct_reorder(Description, -NES_abs),
         label = factor(label, levels = rev(levels(label))))
gsea_data_top_neg <- filter(gsea_data_top, !NES_pos) %>% ungroup %>%
  mutate(label = fct_reorder(Description, -NES_abs))#,
#label = factor(label, levels = rev(levels(label))))

# Positively enriched
xrange_top_pos <- gsea_data_top_pos %>% pull(NES) %>% range %>%
  (function(x) x + top_buffer_x * c(-1,1))
g_top_pos <- ggplot(gsea_data_top_pos, aes(x=NES, y=label, colour=immune_lab)) +
  geom_point() +
  facet_grid(. ~ diversity_order, 
             labeller = labeller(diversity_order = as_labeller(function(x) paste0("q = ", x)))) +
  scale_x_continuous(name = "Normalised Enrichment Score",
                     breaks = seq(-10,10,0.2), limits = xrange_top_pos) +
  scale_colour_immune() +
  theme_stack

# Negatively enriched
xrange_top_neg <- gsea_data_top_neg %>% pull(NES) %>% range %>%
  (function(x) x + top_buffer_x * c(-1,1) * 2)
g_top_neg <- ggplot(gsea_data_top_neg, aes(x=NES, y=label, colour=immune_lab)) +
  geom_point() +
  facet_grid(. ~ diversity_order, 
             labeller = labeller(diversity_order = as_labeller(function(x) paste0("q = ", x)))) +
  scale_x_continuous(name = "Normalised Enrichment Score",
                     breaks = seq(-10,10,0.2), limits = xrange_top_neg) +
  scale_colour_immune() +
  theme_stack

#==============================================================================
# SAVE OUTPUT
#==============================================================================

# Table of robust terms
write_csv(gsea_data_robust_summ, out_path_robust)

# Plots
save_fig <- function(path, plot, plot_height, plot_width = 11, device="png"){
  ggsave(filename=path, plot = plot,
         device = device, width = plot_width,
         height = plot_height, units = "cm", dpi = 320, limitsize=FALSE)
}

save_fig(out_path_top_pos, g_top_pos,
         plot_width * plot_ratio_top, plot_width)
save_fig(out_path_top_neg, g_top_neg,
         plot_width * plot_ratio_top, plot_width)
save_fig(out_path_all_pos, g_all_pos,
         plot_width * plot_ratio_all, plot_width)
save_fig(out_path_all_neg, g_all_neg,
         plot_width * plot_ratio_all, plot_width)
