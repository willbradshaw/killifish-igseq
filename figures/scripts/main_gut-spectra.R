###############################################################################
## MAIN TEXT FIGURE                                                          ##
## Repertoire composition and diversity in intestinal killifish samples      ##
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

# Set parameters
palette <- gut_age_palette
groups  <- gut_age_groups
plot_width <- snakemake@params[["plot_width"]]
plot_height_main <- plot_width * snakemake@params[["plot_ratio_main"]] 
plot_height_pcoa <- plot_width * snakemake@params[["plot_ratio_pcoa"]] 
legend_keysize <- 0.7

# Specify input paths
design_path            <- snakemake@input[["design"]]
clone_alpha_path       <- snakemake@input[["clone_alpha"]]
clone_kw_path          <- snakemake@input[["clone_kw"]]
vj_alpha_path          <- snakemake@input[["vj_alpha"]]
vj_kw_path             <- snakemake@input[["vj_kw"]]
vj_beta_path           <- snakemake@input[["vj_beta"]]
ctab_path              <- snakemake@input[["changeo_gut"]]

#==============================================================================
# Auxiliary functions
#==============================================================================

pull_rdis <- function(dtab, age_group) dtab %>% 
  filter(AGE_WEEKS == age_group) %>% pull(RDI)

get_mwu_grid <- function(dtab){
  age_groups <- unique(dtab[["AGE_WEEKS"]])
  L <- lapply(1:(length(age_groups)-1), function(a) 
    lapply((a+1):length(age_groups), function(b) wilcox.test(
      pull_rdis(dtab, age_groups[a]), pull_rdis(dtab, age_groups[b]))$p.value))
  mwu_grid <- melt(L, value.name = "P") %>%
    group_by(L1) %>% 
    mutate(AGE1 = age_groups[L1], AGE2 = age_groups[L1 + L2]) %>%
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
                               palette = gut_age_palette){
  g <- ggplot(spectrum_data, aes(x=Q)) +
    geom_ribbon(aes(ymin=D_LOWER_EMPIRICAL, ymax=D_UPPER_EMPIRICAL, fill=AGE_WEEKS),
                alpha = ribbon_alpha) +
    geom_line(aes(colour=AGE_WEEKS, y=D), size = line_size) +
    scale_colour_manual(values = palette, name = "Age group (weeks)") +
    scale_fill_manual(values = palette, name = "Age group (weeks)") +
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

#==============================================================================
# A - DESIGN
#==============================================================================

design <- readPNG(design_path)
design_grob <- rasterGrob(design, interpolate=TRUE)
g_design <- ggplot() +
  annotation_custom(design_grob) +
  theme_nothing() + theme(plot.margin = margin(t = 0.1, l = 0.2, r = 0.2, unit="cm"))

fig_a <- g_design

#==============================================================================
# B - CLONAL DIVERSITY
#==============================================================================

# Import data
clone_alpha <- import_div(clone_alpha_path) %>%
  mutate(AGE_WEEKS = factor(AGE_WEEKS, levels = groups)) %>%
  filter(DIVTYPE == "alpha")
clone_kw <- import_div(clone_kw_path) %>% filter(METHOD == "KW permutation test") %>%
  mutate(LABEL = signif_stars(P))

# Make label DB
clone_alpha_labels <- make_label_tab(clone_alpha, clone_kw, 
                                     c(32, 40, 60, 60, 40, 40), 
                                     c(0.5, rep(0.5, 5)))

# Make spectrum plot
g_clone_alpha <- make_spectrum_plot(clone_alpha, clone_alpha_labels,
                                    "Clonal alpha diversity", 500, 500, 100,
                                    1.5)
fig_b <- g_clone_alpha

#==============================================================================
# C - V/J DIVERSITY (ALL)
#==============================================================================

# Import data
vj_alpha <- import_div(vj_alpha_path) %>%
  mutate(AGE_WEEKS = factor(AGE_WEEKS, levels = groups)) %>% 
  filter(DIVTYPE == "alpha")
vj_kw <- import_div(vj_kw_path) %>% filter(METHOD == "KW permutation test") %>%
  mutate(LABEL = signif_stars(P))

# Make label DB
vj_alpha_labels <- make_label_tab(vj_alpha, vj_kw, 
                                  c(7,9,8,8,8,8), 
                                  c(0.3, rep(0.5, 5)))

# Make spectrum plot
g_vj_alpha <- make_spectrum_plot(vj_alpha, vj_alpha_labels,
                                    "VJ alpha diversity", 120, 120, 30,
                                    1.5)
fig_c <- g_vj_alpha

#==============================================================================
# D - BETA SPECTRA
#==============================================================================

# Filter to beta diversity and rescale by sample size
vj_beta <- import_div(vj_beta_path) %>%
  mutate(AGE_WEEKS = factor(AGE_WEEKS, levels = groups)) %>% 
  filter(DIVTYPE == "beta") %>%
  mutate(D = (D-1)/(N_GROUP-1), D_UPPER_EMPIRICAL = (D_UPPER_EMPIRICAL-1)/(N_GROUP-1),
         D_LOWER_EMPIRICAL = (D_LOWER_EMPIRICAL-1)/(N_GROUP-1), D_SD = D_SD/(N_GROUP-1))

# Make spectrum plot
g_vj_beta <- make_spectrum_plot(vj_beta, NULL,
                                "Relative VJ beta diversity", 0.35, 0.3, 0.1, 
                                coord_ratio = 1.5)

fig_d <- g_vj_beta

#==============================================================================
# H - RDI BOXPLOTS
#==============================================================================

cols_gut <- cols(REPLICATE = "c", INDIVIDUAL = "c", FISH = "c")
if (!exists("ctab_gut_age")){
  ctab_gut_age <- suppressMessages(read_tsv(ctab_path, col_types = cols_gut))
}

# Filter ambiguous segment calls
ctab_ambig <- filter(ctab_gut_age, HAS_VJ, !VJ_AMBIG)

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
indivs <- ctab_gut_age %>% group_by(AGE_WEEKS, INDIVIDUAL) %>% summarise %>%
  ungroup() %>% arrange(INDIVIDUAL) %>% mutate(N = row_number())
indiv_index <- setNames(indivs$N, indivs$INDIVIDUAL)
indiv_ages  <-  setNames(indivs$AGE_WEEKS, indivs$INDIVIDUAL)
dist_tab <- rdi %>% as.matrix %>% 
  melt(varnames = c("ID1", "ID2"), value.name = "RDI") %>% as_tibble %>%
  mutate(ROW = indiv_index[ID1], COL = indiv_index[ID2],
         AGE1 = indiv_ages[ID1], AGE2 = indiv_ages[ID2])

# Remove redundant distances and restrict to matching age groups
dist_tab_filtered <- filter(dist_tab, ROW > COL, AGE1 == AGE2) %>%
  select(AGE_WEEKS = AGE1, ID1, ID2, RDI)

# Define maximum distances for annotation positioning
dist_max_tab <- dist_tab_filtered %>% group_by(AGE_WEEKS) %>%
  summarise(DIST_MAX = max(RDI))

# Prepare MWU significance tables
signif_level <- 0.05
mwu_grid_all <- dist_tab_filtered %>% get_mwu_grid %>%
  filter(P <= signif_level) %>%
  inner_join(dist_max_tab %>% rename(AGE1 = AGE_WEEKS), by="AGE1") %>%
  rename(DIST_MAX_AGE1 = DIST_MAX) %>%
  inner_join(dist_max_tab %>% rename(AGE2 = AGE_WEEKS), by="AGE2") %>%
  rename(DIST_MAX_AGE2 = DIST_MAX) %>%
  mutate(DIST_MAX = pmax(DIST_MAX_AGE1,DIST_MAX_AGE2)) %>%
  mutate(DIST_START = min(DIST_MAX)) %>%
  select(-DIST_MAX_AGE1, -DIST_MAX_AGE2) %>%
  mutate(AGE1 = as.numeric(AGE1), AGE2 = as.numeric(AGE2),
         AGE_AVG = (AGE1 + AGE2)/2, AGE_MIN = pmin(AGE1,AGE2),
         AGE_DIFF = pmax(AGE1, AGE2) - pmin(AGE1,AGE2)) %>%
  arrange(AGE_MIN, AGE_DIFF) %>%
  mutate(Y_BAR = DIST_START + 0.5,
         Y_LAB = Y_BAR + 0.2,
         LABEL = signif_stars(P))


# Make RDI boxplot
g_rdi_all <- ggplot(dist_tab_filtered, aes(x=as.numeric(AGE_WEEKS), y=RDI)) +
  geom_boxplot(aes(fill = factor(AGE_WEEKS, levels = groups)), 
               outlier.shape = NA) +
  geom_point(size = 0.7, alpha = 0.2, shape=16) +
  scale_fill_manual(values = gut_age_palette, name = "Age group (days)") +
  scale_colour_manual(values = gut_age_palette, name = "Age group (days)") +
  scale_y_continuous(name="Pairwise RDI distance", limits=c(7, 15),
                     breaks = seq(8, 14, 2)) +
  scale_x_continuous(name="Age group (weeks)", limits=c(3, 19),
                     breaks=c(6, 11, 16)) + 
  coord_fixed((16-6)/(15-7) / 1.5) +
  theme_base + theme(
    legend.position = "none"
    )

g_rdi_all_annot <- g_rdi_all +
  geom_segment(aes(x=AGE1,xend=AGE2,y=Y_BAR,yend=Y_BAR), data = mwu_grid_all,
               size = 0.4) + 
  geom_text(aes(x=AGE_AVG, y=Y_LAB, label=LABEL), data = mwu_grid_all, 
            size=fontsize_base * 5/14 * 1.5)

# Remove boxplot borders
g_rdi_grob <- ggplotGrob(g_rdi_all_annot)
for (n in 1:length(g_rdi_grob$grobs[[6]]$children[[3]]$children)){
  h <- which(grepl("geom_crossbar",names(g_rdi_grob$grobs[[6]]$children[[3]]$children[[n]]$children)))
  g_rdi_grob$grobs[[6]]$children[[3]]$children[[n]]$children[[h]]$children[[1]]$gp$col <- g_rdi_grob$grobs[[6]]$children[[3]]$children[[n]]$children[[h]]$children[[1]]$gp$fill
}

fig_e <- g_rdi_grob

#==============================================================================
# EXTRA - PCOA
#==============================================================================

indivs <- ctab_gut_age %>% group_by(AGE_WEEKS, INDIVIDUAL) %>%
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
pcoa_scale <- 3.2
g_pcoa <- ggplot(pc_tab, aes(x=PCO1, y=PCO2,
                             colour = factor(AGE_WEEKS, levels = groups))) + 
  geom_point(size = 1.5, alpha = 0.8, shape = 16) +
  scale_colour_manual(values = gut_age_palette, name = "Age group (days)") +
  facet_wrap(.~factor(paste(AGE_WEEKS, "weeks"), levels = c("6 weeks", "16 weeks")),
             scales = "fixed") +
  coord_fixed(ratio = 1) +
  scale_x_continuous(name=paste0("Principle co-ord. 1 (", pc_var[1], "%)"),
                     limits = c(-6, 6), breaks = seq(-6, 6, 3)) +
  scale_y_continuous(name=paste0("Principle co-ord. 2 (", pc_var[2], "%)"),
                     limits = c(-4, 8), breaks = seq(-4, 8, 2)) +
  theme_base +
  theme(legend.position = "none",
        panel.grid = element_line(colour = "#EFEFEF"))

#==============================================================================
# LEGEND
#==============================================================================

keysize <- 0.7

g_legend <- clone_alpha %>% filter(Q %in% c(0,1,2)) %>%
  ggplot(aes(x=Q)) +
  geom_ribbon(aes(ymin=D_LOWER_EMPIRICAL, ymax=D_UPPER_EMPIRICAL, fill=AGE_WEEKS),
              alpha = 0.4) +
  geom_line(aes(colour=AGE_WEEKS, y=D), size = 0.5) +
  geom_point(aes(colour=AGE_WEEKS, y=D), size = 1.3) +
  scale_colour_manual(values = gut_age_palette, name = "Age group (weeks)") +
  scale_fill_manual(values = gut_age_palette, name = "Age group (weeks)") +
  guides(colour = guide_legend(keyheight=keysize, keywidth = keysize),
         fill   = guide_legend(keyheight=keysize, keywidth = keysize)) +
  theme_base + theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(margin = margin(r=4))
    )
legend_out <- get_legend(g_legend)

#==============================================================================
# MAKE FIGURE GRID
#==============================================================================

grid_fig <- function(..., ncol = 1, nrow = 1, labels = "AUTO", 
                     axis = "l", vjust = 0.9){
  plot_grid(..., nrow = nrow, ncol = ncol, align = "hv", axis = axis,
            label_size = fontsize_base * fontscale_label,
            label_fontfamily = titlefont, label_colour = "black",
            label_fontface = "plain", labels = labels, vjust=vjust)
}

fig_row1 <- grid_fig(fig_a, ncol = 1, labels="A", vjust = 1)
fig_row2 <- grid_fig(fig_b, fig_c, ncol = 2, labels = c("B", "C"))
fig_row3 <- grid_fig(fig_d, fig_e, ncol = 2, labels = c("D", "E"))
fig <- grid_fig(fig_row1, fig_row2, fig_row3, nrow = 3, 
                 labels = NULL, rel_heights = c(0.75,1,1))
fig_legend <- plot_grid(fig, legend_out, nrow = 2, rel_heights = c(1, 0.03))

#==============================================================================
# SAVE OUTPUT
#==============================================================================

save_fig(snakemake@output[["main"]], fig_legend, plot_height_main, plot_width)
save_fig(snakemake@output[["pcoa"]], g_pcoa, plot_height_pcoa, plot_width)
