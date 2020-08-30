#------------------------------------------------------------------------------
# Packages
#------------------------------------------------------------------------------

load_packages <- function(package_list){
  for (p in package_list){
    suppressMessages(suppressWarnings(library(p, character.only = TRUE)))
  }
}

packages_default <- c("alakazam", "tidyverse", "cowplot", "ape", "png",
                      "grid", "reshape2")

load_packages(packages_default)

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

#------------------------------------------------------------------------------
# Palettes & groups
#------------------------------------------------------------------------------

# Palettes
gut_age_palette <- c("#e78ac3", "#F5C800")
ageing_palette  <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
pilot_palette   <- c("#999999", "#E69F00", "#56B4E9", "#009E73")

# Groups
ageing_groups <- c("39", "56", "73", "128")
gut_age_groups <- c("6", "16")
pilot_groups   <- paste0("2-", seq(3,6))

#------------------------------------------------------------------------------
# Basic functions
#------------------------------------------------------------------------------

import_div <- function(path, col_types = NULL){
  suppressMessages(read_tsv(path, col_types = col_types))
}

signif_stars <- function(p){
  l <- rep("n.s.", length(p))
  l <- ifelse(p <= 0.05, "*", l)
  l <- ifelse(p <= 0.01, "**", l)
  l <- ifelse(p <= 0.001, "***", l)
  return(l)
}

save_fig <- function(path, plot, plot_height, plot_width = 11, device="png"){
  for (d in device){
    ggsave(filename=path, plot = plot,
           device = d, width = plot_width,
           height = plot_height, units = "cm", dpi = 320, limitsize=FALSE)
  }
}

grid_fig <- function(..., ncol = 1, nrow = 1, labels = "AUTO", 
                     axis = "l", vjust = 0.9){
  plot_grid(..., nrow = nrow, ncol = ncol, align = "hv", axis = axis,
            label_size = fontsize_base * fontscale_label,
            label_fontfamily = titlefont, label_colour = "black",
            label_fontface = "plain", labels = labels, vjust=vjust)
}

#------------------------------------------------------------------------------
# IGOR visualisation
#------------------------------------------------------------------------------


# Define events
entropy_events <- tibble(event = c("productive_sequence","nucleotide_sequence", "model", "v_choice", "d_gene", "j_choice", 
                                   "vd_dinucl", "vd_ins", "dj_dinucl", "dj_ins",
                                   "v_3_del", "d_5_del", "d_3_del", "j_5_del"),
                         category = c("Productive sequence","Nucleotide sequence", "Recombination events (total)", rep("Gene choice", 3), 
                                      rep("Insertions",4),rep("Deletions",4)),
                         short = c("productive sequence","nucleotide sequence", "total", "V", "D", "J", 
                                   "VD nts", "VD\nlen", "DJ nts", "DJ\nlen",
                                   "delV", "delD", "delD", "delJ")
)

# Define colour palette
igor_entropy_palette <- c("#FFA500","#FF0000", "#FFCCCC","#66B266","#8080E6","#E66666")
palette_indels <- c("#999999", "#E69F00", "#56B4E9", "#009E73")

import_entropies <- function(path){
  col_entropies <- cols(h = "d", .default = "c")
  htab <- import_div(path, col_entropies) %>% 
    full_join(entropy_events, by = "event") %>%
    group_by(category, short, id) %>% summarise(h = sum(h))
  return(htab)
}

entropy_table <- function(entropies, events = entropy_events){
  # Get levels and references from events table
  categories <- events %>% pull(category) %>% unique
  shorts <- events %>% pull(short) %>% unique %>% 
    (function(x) c(categories, x[-1])) %>% unique
  productive_short <- events %>% filter(event == "productive_sequence") %>% pull(short)
  productive_cat <- events %>% filter(event == "productive_sequence") %>% pull(category)
  model_short <- events %>% filter(event == "model") %>% pull(short)
  model_cat <- events %>% filter(event == "model") %>% pull(category)
  nucleotide_short <- events %>% filter(event == "nucleotide_sequence") %>% pull(short)
  nucleotide_cat <- events %>% filter(event == "nucleotide_sequence") %>% pull(category)
  h_total <- entropies %>% filter(short == "total") %>% pull(h)
  # Construct table
  htab <- entropies %>% group_by(category) %>% 
    summarise(h = sum(h)) %>% mutate(short = category) %>% 
    bind_rows(entropies %>% filter(short != model_short & short != nucleotide_short & short != productive_short)) %>%
    mutate(level = ifelse(category == productive_cat, 1,
                          ifelse(category == nucleotide_cat, 2,
                                 ifelse(category == model_cat, 3,
                                        ifelse(short == category, 4, 5))))) %>%
    mutate(label = ifelse(level == 5, short,
                          paste0(short, ": ", round(h, 0), " bits")),
           short = factor(short, levels = shorts),
           category = factor(category, levels = categories)) %>%
    arrange(short) %>% group_by(level) %>%
    mutate(ch = cumsum(h), chl = lag(ch, default = 0), chx = (ch + chl)/2, 
           chx = ifelse(category %in% c(productive_cat, nucleotide_cat), chx,
                        h_total - chx)) %>%
    group_by(category) %>% select(-ch, -chl) %>%
    mutate(alpha = ifelse(level == 1,1, ifelse(level == 2, 0.8, ifelse(level == 3, 0.7, c(0.4,1)))))
  return(htab)
}

plot_entropy_table <- function(htab, pal = igor_entropy_palette){
  # Create stacked entropy composition plot from prepared entropy table
  g <- ggplot(htab) + theme_minimal() + theme_base + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    plot.margin = margin(t=0,b=0,l=0,r=0),
    panel.border = element_blank(), panel.background = element_blank(), 
    legend.position = "none"
  )
  # Create basic plot structure
  g <- g + geom_col(aes(x=level, y=h, fill=category, alpha=short), width = 1) +
    geom_text(aes(x=level, y = chx, label = label),
              hjust = 0.5, vjust = 0.5, family = font,
              size = fontsize_base * 5/14, lineheight = 0.8) + 
    coord_flip() + scale_x_reverse() + scale_y_reverse()
  # Configure colour display
  g <- g + scale_alpha_manual(values = htab$alpha) +
    scale_fill_manual(values = pal)
  return(g)
}

import_segments <- function(path){
  tab <- import_div(path, col = cols(p = "d", .default = "c")) %>%
    mutate(segtypes = paste0(ifelse(is.na(v), "", "V"),
                             ifelse(is.na(d), "", "D"),
                             ifelse(is.na(j), "", "J"))) %>%
    mutate(v = gsub("_", " / ", v),
           d = gsub("_", " / ", d),
           j = gsub("_", " / ", j))
}
