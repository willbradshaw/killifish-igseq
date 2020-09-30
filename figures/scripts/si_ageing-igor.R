###############################################################################
## SUPPLEMENTARY FIGURE                                                      ##
## Ageing IGoR plots                                                        ##
###############################################################################

#==============================================================================
# Preamble
#==============================================================================
# 
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Source formatting
source("scripts/aux_format-plots.R")

# Specify input paths
v_names_path   <- snakemake@input[["vh"]]
d_names_path   <- snakemake@input[["dh"]]
j_names_path   <- snakemake@input[["jh"]]
segments_path  <- snakemake@input[["segments"]]
indels_path    <- snakemake@input[["indels"]]
entropies_path <- snakemake@input[["entropies"]]
# v_names_path   <- "input_files/data/nfu_vh_name_conversion.csv"
# d_names_path   <- "input_files/data/nfu_dh_name_conversion.csv"
# j_names_path   <- "input_files/data/nfu_jh_name_conversion.csv"
# segments_path  <- "../data_processed/igor/ageing-individual-segments.tsv.gz"
# indels_path    <- "../data_processed/igor/ageing-individual-indels.tsv.gz"
# entropies_path <- "../data_processed/igor/ageing-individual-entropies.tsv.gz"

# Specify run parameters
age_groups  <- ageing_groups
age_palette <- ageing_palette
trace_alpha <- 0.7
trace_size  <- 0.4
trace_linetype <- "solid"
keysize_legend <- 0.7

# Specify output parameters
out_path_indels <- snakemake@output[["indels"]]
out_path_segments <- snakemake@output[["segments"]]
out_path_entropies <- snakemake@output[["entropies"]]
out_width <- snakemake@params[["plot_width"]]
out_height_indels <- snakemake@params[["plot_ratio_indels"]] * out_width
out_height_entropies <- snakemake@params[["plot_ratio_entropies"]] * out_width
out_height_segments <- snakemake@params[["plot_ratio_segments"]] * out_width
# out_path_indels <- "output_files/dev_si_ageing-igor-indels.png"
# out_path_segments <- "output_files/dev_si_ageing-igor-segments.png"
# out_path_entropies <- "output_files/dev_si_ageing-igor-entropies.png"
# out_width <- 11
# out_height_indels <- out_width * 0.8
# out_height_segments <- out_width * 0.8
# out_height_entropies <- out_width * 0.6

# Modify theme
theme_base <-   theme_base + theme(
  legend.title = element_text(margin = margin(r = 0.5, unit = "cm")),
  legend.text = element_text(margin = margin(l = 0.1, r = 0.25, unit = "cm")),
  legend.margin = margin(t=0.05, b = 0.15, unit = "cm"),
  # plot.margin = margin(b=0.1, unit = "cm")
)

cat("done.\n")

#==============================================================================
# Import data
#==============================================================================

cat("\nPreparing data...")

# Segment name correspondence
v_names <- suppressMessages(read_csv(v_names_path))
d_names <- suppressMessages(read_csv(d_names_path))
j_names <- suppressMessages(read_csv(j_names_path))

# Segments
col_segments <- cols(p = "d", .default = "c")
segments <- import_div(segments_path, col = col_segments) %>%
  mutate(SEGMENT_TYPES = paste0(ifelse(is.na(v), "", "V"),
                                ifelse(is.na(d), "", "D"),
                                ifelse(is.na(j), "", "J"))) %>%
  rename(ID = id, P = p, V_NAME_NEW = v, D_NAME_NEW = d, J_NAME_NEW = j) %>%
  full_join(v_names, by = "V_NAME_NEW") %>%
  full_join(d_names, by = "D_NAME_NEW") %>%
  full_join(j_names, by = "J_NAME_NEW") %>%
  mutate(V_NAME_OLD = gsub("_", " / ", V_NAME_OLD),
         D_NAME_OLD = gsub("_", " / ", D_NAME_OLD),
         J_NAME_OLD = gsub("_", " / ", J_NAME_OLD)) %>%
  mutate(AGE_GROUP = as.integer(sub("-.*", "", ID)),
         AGE_DAYS = age_groups[AGE_GROUP])

# Indels
col_indels <- cols(p = "d", n = "i", .default = "c")
indels <- import_div(indels_path, col_indels) %>%
  mutate(age_group = as.integer(sub("-.*", "", id)),
         age_days = age_groups[age_group])

cat("done.\n")

#==============================================================================
# Make segment plots
#==============================================================================

cat("\nPreparing segment plots...")

# Define function
plot_segments <- function(segment_type, keysize = keysize_legend){
  g <- segments %>% filter(SEGMENT_TYPES == segment_type) %>%
    ggplot(aes_string(x=paste0(segment_type, "_NAME_OLD"), y="P",
                      colour = "AGE_DAYS", group = "ID")) +
    geom_line(alpha = trace_alpha, size = trace_size,
              linetype = trace_linetype) +
    scale_colour_manual(values = age_palette, name = "Age group (days)") +
    scale_y_continuous(name= "Probability (%)", labels = function(y) y * 100,
                       expand = c(0,0)) +
    guides(colour = guide_legend(keyheight=keysize, keywidth = keysize),
           fill   = guide_legend(keyheight=keysize, keywidth = keysize)) +
    theme_base + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = fontsize_base * fontscale_title,
                                 angle = 45, hjust = 1),
      legend.title = element_text(margin = margin(r = 0.5, unit = "cm")),
      legend.text = element_text(margin = margin(l = 0.1, r = 0.25, unit = "cm")),
    )
  return(g)
}

# Define individual plots
g_segment_v <- plot_segments("V")
g_segment_d <- plot_segments("D")
g_segment_j <- plot_segments("J")

# Extract legend
legend_segment <- get_legend(g_segment_v)

# Prepare grid
g_segment_dj <- plot_grid(g_segment_d + theme(legend.position = "none"), 
                          g_segment_j + theme(legend.position = "none"), 
                          nrow = 1, ncol = 2,
                          align = "hv", axis = "l", labels = c("B", "C"),
                          label_size = fontsize_base * fontscale_label,
                          label_fontfamily = titlefont, label_colour = "black",
                          label_fontface = "plain", vjust=1.1)
g_segment_vv <- plot_grid(g_segment_v + theme(legend.position = "none"),
                          nrow = 1, ncol = 1,
                          align = "hv", axis = "l", labels = c("A"),
                          label_size = fontsize_base * fontscale_label,
                          label_fontfamily = titlefont, label_colour = "black",
                          label_fontface = "plain", vjust=1.1)
g_segments <- plot_grid(g_segment_vv, g_segment_dj, legend_segment,
                       nrow = 3, ncol = 1,
                       rel_heights = c(1,1,0.08))

cat("done.\n")

#==============================================================================
# Make indel plots
#==============================================================================

cat("\nPreparing indel plots...")

# Prepare event mapping
event_mapping <- list(d_3_del = "3' D-deletions",
                      d_5_del = "5' D-deletions",
                      dj_ins  = "DJ-insertions",
                      j_5_del = "5' J-deletions",
                      v_3_del = "3' V-deletions",
                      vd_ins  = "VD-insertions")
event_order <- c("VD-insertions", "3' V-deletions", "5' D-deletions",
                 "DJ-insertions", "3' D-deletions", "5' J-deletions")
indels_map <- indels %>% 
  mutate(event_long = unlist(event_mapping[event])) %>%
  mutate(event_long = factor(event_long, levels = event_order))

# Make plot
keysize = keysize_legend
g_indels <- ggplot(indels_map, aes(x=n, y=p, colour = age_days,
                                   group = id)) +
  geom_line(alpha = trace_alpha, size = trace_size,
            linetype = trace_linetype) +
  facet_wrap(~event_long, nrow = 2, scales = "free") +
  scale_colour_manual(values = age_palette, name = "Age group (days)") +
  scale_y_continuous(name= "Probability (%)", labels = function(y) y * 100,
                     expand = c(0,0)) +
  scale_x_continuous(name = "# Insertions/Deletions") +
  guides(colour = guide_legend(keyheight=keysize, keywidth = keysize),
         fill   = guide_legend(keyheight=keysize, keywidth = keysize)) +
  theme_base

cat("done.\n")

#==============================================================================
# Make entropy plots
#==============================================================================

cat("\nPreparing entropy plots...")

# Import entropy data
h_indiv <- import_entropies(entropies_path) %>%
  mutate(age_group = as.integer(sub("-.*", "", id)),
         age_days = as.integer(age_groups[age_group]))

# Aggregate by individual
h_parts <- h_indiv %>% filter(category %in% c("Insertions", "Deletions",
                                              "Gene choice")) %>%
  group_by(id, age_days, category) %>% summarise(h = sum(h))

# Test for age effects
h_pvalues <- h_parts %>% group_by(category) %>%
  summarise(p = kruskal.test(formula = h~age_days)$p.value)

# Make plot
g_entropies <- ggplot(h_parts) +
  geom_boxplot(aes(x=age_days, y=h, fill = factor(age_days, levels = age_groups),
                   group = age_days), outlier.shape = NA) +
  geom_point(aes(x=age_days, y=h), size = 0.5, shape = 16, alpha = 0.5) +
  geom_text(aes(label=paste0("p = ", round(p, 2))), 
            data = h_pvalues,
            x=128, y=12.2, family = font, size = fontsize_base * 5/14,
            hjust=1, vjust=0) +
  scale_fill_manual(values = age_palette, name = "Age (days)") +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Generative entropy (bits)") +
  facet_wrap(~factor(category, levels = c("Gene choice", "Insertions",
                                          "Deletions")), scales = "fixed") +
  scale_colour_manual(values = age_palette, name = "Age group (days)") +
  guides(colour = guide_legend(keyheight=keysize, keywidth = keysize),
         fill   = guide_legend(keyheight=keysize, keywidth = keysize)) +
  theme_base + theme(
    panel.grid = element_line(colour = "#F7F7F7")
  )

# Remove boxplot borders
g_entropies_grob <- ggplotGrob(g_entropies)
gtrees <- which(sapply(g_entropies_grob$grobs, function(g) "gTree" %in% class(g)))

for (n in 1:length(g_entropies_grob$grobs)){
  if (!("gTree" %in% class(g_entropies_grob$grobs[[n]]))) next
  if (length(g_entropies_grob$grobs[[n]]$children) == 0) next
  children_boxplot <- grep("boxplot", names(g_entropies_grob$grobs[[n]]$children))
  if (sum(children_boxplot) == 0) next
  for (m in children_boxplot){
    grandchildren_boxplot <- grep("boxplot", names(g_entropies_grob$grobs[[n]]$children[[m]]$children))
    if (sum(grandchildren_boxplot) == 0) next
    for (p in grandchildren_boxplot){
      h <- grep("geom_crossbar",names(g_entropies_grob$grobs[[n]]$children[[m]]$children[[p]]$children))
      g_entropies_grob$grobs[[n]]$children[[m]]$children[[p]]$children[[h]]$children[[1]]$gp$col <- g_entropies_grob$grobs[[n]]$children[[m]]$children[[p]]$children[[h]]$children[[1]]$gp$fill
    }
  }
}

cat("done.\n")

#==============================================================================
# SAVE OUTPUT
#==============================================================================

cat("\nSaving output...")
save_fig(out_path_indels, g_indels, out_height_indels, out_width)
save_fig(out_path_segments, g_segments, out_height_segments, out_width)
save_fig(out_path_entropies, g_entropies_grob, out_height_entropies, out_width)
cat("...done.\n")