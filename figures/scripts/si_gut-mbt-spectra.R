###############################################################################
## SUPPLEMENTARY FIGURE                                                      ##
## Microbiota-transfer diversity spectra                                     ##
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

# Install RDI from minicran
repo <- normalizePath(snakemake@input[["repo_dir"]]) %>% paste0("file://", .)
lib <- snakemake@input[["lib_dir"]]
suppressMessages(suppressWarnings(install.packages("rdi", repos = repo, 
                                                   lib = lib)))
library(rdi, lib.loc = lib)

# Define group mappingls
mapping <- list(ABX_16 = "16 weeks,\nantibiotic treated",
                SMT_16 = "16 weeks,\nsame-age transfer",
                YMT_16 = "16 weeks,\nyoung transfer",
                YI_6 = "6 weeks,\nuntreated",
                WT_16 = "16 weeks,\nuntreated",
                A_16 = "All 16-week\ngroups"
)

# Define group order and colours
groups_all <- c("WT_16", "ABX_16", "SMT_16", "YMT_16")
groups_long_all <- unlist(mapping[groups_all])
group_palette <- c("#8da0cb", "#e5c494", "#fc8d62", "#a6d854")

# Configure input paths
clone_alpha_path       <- snakemake@input[["clone_alpha"]]
clone_kw_path          <- snakemake@input[["clone_kw"]]
vj_alpha_path          <- snakemake@input[["vj_alpha"]]
vj_kw_path             <- snakemake@input[["vj_kw"]]
vj_beta_path           <- snakemake@input[["vj_beta"]]
ctab_path              <- snakemake@input[["changeo_db"]]

# Define output parameters
out_path <- snakemake@output[[1]]
out_width <- snakemake@params[["plot_width"]]
out_height <- snakemake@params[["plot_ratio"]] * out_width

# Modify theme
theme_base <- theme_base + theme(plot.margin = margin(b = 0, t = 0,
                                                      l = 0.2, r = 0.2, unit = "cm"))

#==============================================================================
# Auxiliary functions
#==============================================================================

pull_rdis <- function(dtab, age_group) dtab %>% 
  filter(GROUP_LONG == age_group) %>% pull(RDI)

get_mwu_grid <- function(dtab){
  age_groups <- unique(dtab[["GROUP_LONG"]])
  L <- lapply(1:(length(age_groups)-1), function(a) 
    lapply((a+1):length(age_groups), function(b) wilcox.test(
      pull_rdis(dtab, age_groups[a]), pull_rdis(dtab, age_groups[b]))$p.value))
  mwu_grid <- melt(L, value.name = "P") %>%
    group_by(L1) %>% 
    mutate(GRP1 = age_groups[L1], GRP2 = age_groups[L1 + L2]) %>%
    ungroup() %>% select(-L1, -L2)
  return(mwu_grid)
}

make_label_tab <- function(spectrum_data, kw_data, ystep, hjust,
                           star_scale = 1.5,
                           qvals = c(0,1,1.5,2,3,4)){
  spectrum_data %>% filter(Q %in% qvals) %>% group_by(Q) %>% 
    summarise(DMIN = min(D), DMAX = max(D)) %>% mutate(DLAB = DMAX + ystep) %>%
    inner_join(kw_data, by="Q") %>%
    mutate(SIZE = ifelse(LABEL == "n.s.", 1, star_scale),
           VJUST = ifelse(LABEL == "n.s.", 0, 0.5),
           HJUST = hjust)
}

make_spectrum_plot <- function(spectrum_data, label_data, ylab,
                               ylim_max, ybreak_max, ybreak_step, 
                               coord_ratio = 1, xlim_max = 4.05, xbreak_max = 4,
                               ribbon_alpha = 0.4, line_size = 0.3,
                               segment_size = 0.4, segment_linetype = "dotted",
                               palette = group_palette){
  g <- ggplot(spectrum_data, aes(x=Q)) +
    geom_ribbon(aes(ymin=D_LOWER_EMPIRICAL, ymax=D_UPPER_EMPIRICAL, fill=GROUP_LONG),
                alpha = ribbon_alpha) +
    geom_line(aes(colour=GROUP_LONG, y=D), size = line_size) +
    scale_colour_manual(values = palette, name = "Age group (days)") +
    scale_fill_manual(values = palette, name = "Age group (days)") +
    scale_x_continuous(name="Diversity order (q)", #expand = c(0,0),
                       limits = c(0, xlim_max), 
                       breaks = seq(0, xbreak_max, 1)) +
    scale_y_continuous(name=ylab, expand = c(0,0),
                       limits = c(0, ylim_max),
                       breaks = seq(0, ybreak_max, ybreak_step)) +
    coord_fixed(ratio = xlim_max/ylim_max / coord_ratio) +
    theme_base + theme(legend.position = "none")
  if (!is.null(label_data)){
    g <- g + geom_segment(aes(y=DMIN, yend=DMAX, xend=Q), 
                 data = label_data, linetype = segment_linetype,
                 size = segment_size) +
      geom_text(aes(x=Q, y=DLAB, label=LABEL), data = label_data,
                size = label_data$SIZE * fontsize_base * 5/14, vjust = label_data$VJUST,
                hjust = label_data$HJUST)
  }
  return(g)
}

cat("done.\n")

#==============================================================================
# A - CLONAL DIVERSITY
#==============================================================================

cat("\nGenerating clonal alpha-diversity plot...")

# Import data
clone_alpha <- import_div(clone_alpha_path) %>%
  mutate(GROUP_LONG = unlist(mapping[GROUP])) %>%
  mutate(GROUP_LONG = factor(GROUP_LONG, levels = groups_long_all)) %>%
  filter(GROUP != "YI_6")
clone_kw <- import_div(clone_kw_path) %>% filter(METHOD == "KW permutation test") %>%
  mutate(LABEL = signif_stars(P))

# Make label DB
clone_alpha_labels <- make_label_tab(clone_alpha, clone_kw, 
                                     c(15, 20, 15, 15, 15, 15), 
                                     c(0.2, rep(0.5, 5)))

# Make spectrum plot
g_clone_alpha <- make_spectrum_plot(clone_alpha, clone_alpha_labels,
                                    "Clonal alpha diversity", 240, 240, 80,
                                    1.5)

cat("done.\n")

#==============================================================================
# B - V/J DIVERSITY
#==============================================================================

cat("\nGenerating VJ alpha-diversity plot...")

# Import data
vj_alpha <- import_div(vj_alpha_path) %>%
  mutate(GROUP_LONG = unlist(mapping[GROUP])) %>%
  mutate(GROUP_LONG = factor(GROUP_LONG, levels = groups_long_all)) %>%
  filter(GROUP != "YI_6")
vj_kw <- import_div(vj_kw_path) %>% filter(METHOD == "KW permutation test") %>%
  mutate(LABEL = signif_stars(P))

# Make label DB
vj_alpha_labels <- make_label_tab(vj_alpha, vj_kw, 
                                     c(9, 9, 6, 5, 4, 4), 
                                     c(0.2, rep(0.5, 5)))

# Make spectrum plot
g_vj_alpha <- make_spectrum_plot(vj_alpha, vj_alpha_labels,
                                    "VJ alpha diversity", 80, 80, 20,
                                    1.5)

cat("done.\n")

#==============================================================================
# C - BETA SPECTRA
#==============================================================================

cat("\nGenerating VJ beta-diversity plot...")

# Filter to beta diversity and rescale by sample size
vj_beta <- import_div(vj_beta_path) %>%
  mutate(GROUP_LONG = unlist(mapping[GROUP])) %>%
  mutate(GROUP_LONG = factor(GROUP_LONG, levels = groups_long_all)) %>%
  filter(GROUP != "YI_6") %>%
  mutate(D = (D-1)/(N_GROUP-1), D_UPPER_EMPIRICAL = (D_UPPER_EMPIRICAL-1)/(N_GROUP-1),
         D_LOWER_EMPIRICAL = (D_LOWER_EMPIRICAL-1)/(N_GROUP-1), D_SD = D_SD/(N_GROUP-1))

# Make spectrum plot
g_vj_beta <- make_spectrum_plot(vj_beta, NULL,
                                "Relative VJ beta diversity", 1.05, 1, 0.2)

cat("done.\n")

#==============================================================================
# D - RDI BOXPLOTS
#==============================================================================

cat("\nGenerating RDI plot...")

if (!exists("ctab_gut_group")){
  ctab_gut_group <- import_div(ctab_path, cols(REPLICATE = "c", 
                                               INDIVIDUAL = "c", FISH = "c")) %>%
    mutate(GROUP_LONG = unlist(mapping[GROUP])) %>%
    mutate(GROUP_LONG = factor(GROUP_LONG, levels = groups_long_all)) %>%
    filter(GROUP != "YI_6")
}

# Filter ambiguous segment calls
ctab_ambig <- filter(ctab_gut_group, HAS_VJ, !VJ_AMBIG)

# Extract gene calls and annotations
genes <- pull(ctab_ambig, BEST_VJ_CALL)
annots <- pull(ctab_ambig, INDIVIDUAL)

# Compute RDI segment calls
counts <- calcVDJcounts(genes = genes, seqAnnot = annots,
                        simplifyNames = TRUE, splitCommas = FALSE)

# Compute RDI distance matrix
rdi <- calcRDI(counts, subsample = TRUE, distMethod = "euclidean", nIter = 100,
               constScale = TRUE, units = "lfc")

# Convert RDI output to distance tibbles
indivs <- ctab_gut_group %>% group_by(GROUP_LONG, INDIVIDUAL) %>% summarise %>%
  ungroup() %>% arrange(INDIVIDUAL) %>% mutate(N = row_number())
indiv_index <- setNames(indivs$N, indivs$INDIVIDUAL)
indiv_groups  <-  setNames(indivs$GROUP_LONG, indivs$INDIVIDUAL)
dist_tab <- rdi %>% as.matrix %>% 
  melt(varnames = c("ID1", "ID2"), value.name = "RDI") %>% as_tibble %>%
  mutate(ROW = indiv_index[ID1], COL = indiv_index[ID2],
         GRP1 = indiv_groups[ID1], GRP2 = indiv_groups[ID2])

# Remove redundant distances and restrict to matching age groups
dist_tab_filtered <- filter(dist_tab, ROW > COL, GRP1 == GRP2) %>%
  select(GROUP_LONG = GRP1, ID1, ID2, RDI)

# Define maximum distances for annotation positioning
dist_max_tab <- dist_tab_filtered %>% group_by(GROUP_LONG) %>%
  summarise(DIST_MAX = max(RDI))

# Prepare MWU significance tables
signif_level <- 0.05
mwu_grid_all <- dist_tab_filtered %>% get_mwu_grid %>%
  filter(P <= signif_level) %>%
  inner_join(dist_max_tab %>% rename(GRP1 = GROUP_LONG), by="GRP1") %>%
  rename(DIST_MAX_GRP1 = DIST_MAX) %>%
  inner_join(dist_max_tab %>% rename(GRP2 = GROUP_LONG), by="GRP2") %>%
  rename(DIST_MAX_GRP2 = DIST_MAX) %>%
  mutate(DIST_MAX = pmax(DIST_MAX_GRP1,DIST_MAX_GRP2)) %>%
  mutate(DIST_START = min(DIST_MAX)) %>%
  select(-DIST_MAX_GRP1, -DIST_MAX_GRP2) %>%
  mutate(GRP1 = as.numeric(GRP1), GRP2 = as.numeric(GRP2),
         AGE_AVG = (GRP1 + GRP2)/2, AGE_MIN = pmin(GRP1,GRP2),
         AGE_DIFF = pmax(GRP1, GRP2) - pmin(GRP1,GRP2)) %>%
  arrange(AGE_MIN, AGE_DIFF) %>%
  mutate(Y_BAR = DIST_START + c(0.2, 0.6, 1.0),
         Y_LAB = Y_BAR + 0.1,
         LABEL = signif_stars(P))

# Make RDI boxplot
g_rdi_all <- ggplot(dist_tab_filtered, aes(x=GROUP_LONG, y=RDI)) +
  geom_boxplot(aes(fill = GROUP_LONG), 
               outlier.shape = NA) +
  geom_point(size = 0.7, alpha = 0.2, shape=16) +
  scale_fill_manual(values = group_palette, name = "Treatment group") +
  scale_colour_manual(values = group_palette, name = "Treatment group") +
  scale_y_continuous(name="Pairwise RDI distance", limits=c(9,14),
                     breaks = seq(9,14,1), expand = c(0,0)) +
  scale_x_discrete(name="Treatment group") + 
  coord_fixed(4/(14-9)) +
  theme_base + theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    )

if (nrow(mwu_grid_all) > 0){
  g_rdi_all_annot <- g_rdi_all +
    geom_segment(aes(x=GRP1,xend=GRP2,y=Y_BAR,yend=Y_BAR), data = mwu_grid_all,
                 size = 0.4) + 
    geom_text(aes(x=AGE_AVG, y=Y_LAB, label=LABEL), data = mwu_grid_all, 
              size=fontsize_base * 5/14 * 1.5)
} else {
  g_rdi_all_annot <- g_rdi_all
}

# Remove boxplot borders
g_rdi_grob <- ggplotGrob(g_rdi_all_annot)
for (n in 1:length(g_rdi_grob$grobs[[6]]$children[[3]]$children)){
  h <- which(grepl("geom_crossbar",names(g_rdi_grob$grobs[[6]]$children[[3]]$children[[n]]$children)))
  g_rdi_grob$grobs[[6]]$children[[3]]$children[[n]]$children[[h]]$children[[1]]$gp$col <- g_rdi_grob$grobs[[6]]$children[[3]]$children[[n]]$children[[h]]$children[[1]]$gp$fill
}

cat("done.\n")

#==============================================================================
# E - PCOA
#==============================================================================

cat("\nGenerating PCoA plot...")

indivs <- ctab_gut_group %>% group_by(GROUP_LONG, INDIVIDUAL) %>% 
  summarise %>% ungroup

# Perform PCoA
pc <- pcoa(rdi)

# Extract co-ordinates in first two axes into tibble for plotting
pc_tab <- tibble(INDIVIDUAL = rownames(as.matrix(rdi)),
                 PCO1 = pc$vectors[,1], PCO2 = pc$vectors[,2]) %>%
  full_join(indivs, by = "INDIVIDUAL")

# Obtain proportion of variance in each dimension
pc_var <- pc$values %>% pull(Broken_stick) %>% (function(x) round(x*100, 1))

# Make PCoA plot
pcoa_scale <- 8

g_pcoa <- ggplot(pc_tab, aes(x=PCO1, y=PCO2,
                               colour = GROUP_LONG)) + 
  geom_point(size = 1, alpha = 0.8, shape = 16) +
  scale_colour_manual(values = group_palette, name = "Treatment group") +
  facet_wrap(~GROUP_LONG, scales = "fixed", ncol = 2) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=axis_width) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=axis_width) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(name=paste0("Principle co-ord. 1 (", pc_var[1], "%)"),
                     limits = c(-pcoa_scale,pcoa_scale),
                     breaks = seq(-pcoa_scale,pcoa_scale,pcoa_scale/2)) +
  scale_y_continuous(name=paste0("Principle co-ord. 2 (", pc_var[2], "%)"),
                     limits = c(-pcoa_scale,pcoa_scale),
                     breaks = seq(-pcoa_scale,pcoa_scale,pcoa_scale/2)) +
  theme_base +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

cat("done.\n")

#==============================================================================
# LEGEND
#==============================================================================

cat("\nPreparing legend...")

keysize <- 0.9

g_legend <- clone_alpha %>% filter(Q %in% c(0,1,2)) %>%
  ggplot(aes(x=Q)) +
  geom_ribbon(aes(ymin=D_LOWER_EMPIRICAL, ymax=D_UPPER_EMPIRICAL, fill=GROUP_LONG),
              alpha = 0.4) +
  geom_line(aes(colour=GROUP_LONG, y=D), size = 0.5) +
  geom_point(aes(colour=GROUP_LONG, y=D), size = 1.3) +
  scale_colour_manual(values = group_palette, name = "Treatment group") +
  scale_fill_manual(values = group_palette, name = "Treatment group") +
  guides(colour = guide_legend(keyheight=keysize, keywidth = keysize),
         fill   = guide_legend(keyheight=keysize, keywidth = keysize)) +
  theme_base + theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(margin = margin(r=4))
    )
legend_out <- get_legend(g_legend)

cat("done.\n")

#==============================================================================
# MAKE FIGURE GRID
#==============================================================================

cat("\nGenerating final grid plot...")

grid_fig <- function(..., ncol = 1, nrow = 1, labels = "AUTO", 
                     axis = "l", vjust = 1.1){
  plot_grid(..., nrow = nrow, ncol = ncol, align = "hv", axis = axis,
            label_size = fontsize_base * fontscale_label,
            label_fontfamily = titlefont, label_colour = "black",
            label_fontface = "plain", labels = labels, vjust=vjust)
}

fig_row1 <- grid_fig(g_clone_alpha, g_vj_alpha, ncol = 2,
                     labels = c("A", "B"))
fig_row2 <- grid_fig(g_vj_beta, g_rdi_grob, g_pcoa, ncol = 3,
                     labels = c("C", "D", "E"), axis = "lb")

fig <- grid_fig(fig_row1, fig_row2, nrow = 2, 
                 labels = NULL, rel_heights = c(1,1))
fig_legend <- plot_grid(fig, legend_out, nrow = 2,
                        rel_heights = c(1, 0.06))

cat("done.\n")

#==============================================================================
# SAVE OUTPUT
#==============================================================================

cat("\nSaving output...")
save_fig(out_path, fig_legend, out_height, out_width)
cat("...done.\n")