###############################################################################
## MAIN TEXT FIGURE                                                          ##
## Repertoire composition and diversity in whole-body killifish samples      ##
###############################################################################

#------------------------------------------------------------------------------
# Preamble
#------------------------------------------------------------------------------

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
age_palette <- ageing_palette
groups  <- ageing_groups
ci_effect <- snakemake@params[["ci_effect"]]
plot_width <- snakemake@params[["plot_width"]]
plot_height <- plot_width * snakemake@params[["plot_ratio"]] 
legend_keysize <- 0.7

# Specify input paths
clone_alpha_path       <- snakemake@input[["clone_alpha"]]
clone_kw_path          <- snakemake@input[["clone_kw"]]
vj_alpha_path          <- snakemake@input[["vj_all_alpha"]]
vj_kw_path             <- snakemake@input[["vj_all_kw"]]
vj_large_alpha_path    <- snakemake@input[["vj_large_alpha"]]
vj_large_kw_path       <- snakemake@input[["vj_large_kw"]]
vj_small_alpha_path    <- snakemake@input[["vj_small_alpha"]]
vj_small_kw_path       <- snakemake@input[["vj_small_kw"]]
vj_small_alpha_bs_path <- snakemake@input[["vj_small_alpha_bs"]]
vj_large_alpha_bs_path <- snakemake@input[["vj_large_alpha_bs"]]
vj_beta_path           <- snakemake@input[["vj_all_beta"]]
ctab_path              <- snakemake@input[["changeo_ageing"]]

#==============================================================================
# Auxiliary functions
#==============================================================================

pull_rdis <- function(dtab, age_group) dtab %>% 
  filter(AGE_DAYS == age_group) %>% pull(RDI)

get_mwu_grid <- function(dtab){
  ages <- unique(dtab[["AGE_DAYS"]])
  L <- lapply(1:(length(ages)-1), function(a) 
    lapply((a+1):length(ages), function(b) wilcox.test(
      pull_rdis(dtab, ages[a]), pull_rdis(dtab, ages[b]))$p.value))
  mwu_grid <- melt(L, value.name = "P") %>%
    group_by(L1) %>% 
    mutate(AGE1 = ages[L1], AGE2 = ages[L1 + L2]) %>%
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
                               palette = age_palette){
  g <- ggplot(spectrum_data, aes(x=Q)) +
    geom_ribbon(aes(ymin=D_LOWER_EMPIRICAL, ymax=D_UPPER_EMPIRICAL, fill=AGE_DAYS),
                alpha = ribbon_alpha) +
    geom_line(aes(colour=AGE_DAYS, y=D), size = line_size) +
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
# A - DESIGN
#==============================================================================

cat("\nGenerating experimental design...")

design_path <- snakemake@input[["design"]]
design <- readPNG(design_path)
design_grob <- rasterGrob(design, interpolate=TRUE)
g_design <- ggplot() +
  annotation_custom(design_grob) +
  theme_nothing() + 
  theme(plot.margin = margin(t = 0.3, l = 0.2, r = 0.2, unit="cm"))

fig_a <- g_design

cat("done.\n")

#==============================================================================
# B - CLONAL DIVERSITY
#==============================================================================

cat("\nGenerating clonal alpha-diversity...")

# Import data
clone_alpha <- import_div(clone_alpha_path) %>%
  mutate(AGE_DAYS = factor(AGE_DAYS, levels = groups)) %>%
  filter(DIVTYPE == "alpha")
clone_kw <- import_div(clone_kw_path) %>% 
  filter(METHOD == "KW permutation test") %>%
  mutate(LABEL = signif_stars(P))

# Make label DB
clone_alpha_labels <- make_label_tab(clone_alpha, clone_kw, 
                                     c(80, 120, 170, 70, 70, 70), 
                                     c(0.2, rep(0.5, 5)))

# Make spectrum plot
g_clone_alpha <- make_spectrum_plot(clone_alpha, clone_alpha_labels,
                                    "Clonal alpha diversity", 1380, 1200, 300,
                                    1.5)
fig_b <- g_clone_alpha

cat("done.\n")

#==============================================================================
# C - V/J DIVERSITY (ALL)
#==============================================================================

cat("\nGenerating VJ alpha-diversity (all clones)...")

# Import data
vj_alpha <- import_div(vj_alpha_path) %>%
  mutate(AGE_DAYS = factor(AGE_DAYS, levels = groups)) %>%
  filter(DIVTYPE == "alpha")
vj_kw <- import_div(vj_kw_path) %>% 
  filter(METHOD == "KW permutation test") %>%
  mutate(LABEL = signif_stars(P))

# Make label DB
vj_alpha_labels <- make_label_tab(vj_alpha, vj_kw, 
                                  c(7, 12, 12, 10, 8, 8), 
                                  c(0.2, rep(0.5, 5)))

# Make spectrum plot
g_vj_alpha <- make_spectrum_plot(vj_alpha, vj_alpha_labels,
                                    "VJ alpha diversity", 120, 120, 30,
                                    1.5)
fig_c <- g_vj_alpha

cat("done.\n")

#==============================================================================
# D - V/J DIVERSITY (LARGE)
#==============================================================================

cat("\nGenerating VJ alpha-diversity (large clones)...")

# Import data
vj_large_alpha <- import_div(vj_large_alpha_path) %>%
  mutate(AGE_DAYS = factor(AGE_DAYS, levels = groups)) %>%
  filter(DIVTYPE == "alpha")
vj_large_kw <- import_div(vj_large_kw_path) %>% 
  filter(METHOD == "KW permutation test") %>%
  mutate(LABEL = signif_stars(P))

# Make label DB
vj_large_alpha_labels <- make_label_tab(vj_large_alpha, vj_large_kw, 
                                  rep(1.5, 6), c(0.5, 0.5, 0.5, 0.2, 0.5, 0.5))

# Make spectrum plot
g_vj_large_alpha <- make_spectrum_plot(vj_large_alpha, vj_large_alpha_labels,
                                 "VJ α-div, large clones", 24, 24, 6)

fig_d <- g_vj_large_alpha

cat("done.\n")

#==============================================================================
# E - V/J DIVERSITY (SMALL)
#==============================================================================

cat("\nGenerating VJ alpha-diversity (small clones)...")

# Import data
vj_small_alpha <- import_div(vj_small_alpha_path) %>%
  mutate(AGE_DAYS = factor(AGE_DAYS, levels = groups)) %>%
  filter(DIVTYPE == "alpha")
vj_small_kw <- import_div(vj_small_kw_path) %>% 
  filter(METHOD == "KW permutation test") %>%
  mutate(LABEL = signif_stars(P))

# Make label DB
vj_small_alpha_labels <- make_label_tab(vj_small_alpha, vj_small_kw, 
                                        c(7, 12, 10, 8, 6, 8), 
                                        c(0.1, 0.5, 0.5, 0.5, 0.5, 0.5))

# Make spectrum plot
g_vj_small_alpha <- make_spectrum_plot(vj_small_alpha, vj_small_alpha_labels,
                                       "VJ α-div, small clones", 120, 120, 30)

fig_e <- g_vj_small_alpha

cat("done.\n")

#==============================================================================
# F - EFFECT SIZES
#==============================================================================

cat("\nGenerating effect-size plot...")

vj_small_alpha_bs <- import_div(vj_small_alpha_bs_path) %>% mutate(SIZE="small")
vj_large_alpha_bs <- import_div(vj_large_alpha_bs_path) %>% mutate(SIZE="large")
vj_alpha_bs <- bind_rows(vj_small_alpha_bs, vj_large_alpha_bs)

# Process
vj_alpha_bs_proc <- vj_alpha_bs %>% select(-N_GROUP) %>%
  spread(AGE_DAYS, D) %>% 
  mutate(`73 vs. 39` = `73` / `39`,
         `128 vs. 39` = `128` / `39`) %>%
  select(ITER, Q, SIZE, `73 vs. 39`, `128 vs. 39`) %>%
  gather(METRIC, D, -ITER, -Q, -SIZE)

# Compute means & CIs
vj_alpha_bs_ci <- vj_alpha_bs_proc %>%
  group_by(Q, SIZE, METRIC) %>%
  summarise(D_MEAN = mean(D),
            D_LOWER_EMPIRICAL = quantile(D, (1-ci_effect)/2),
            D_UPPER_EMPIRICAL = quantile(D, 1-(1-ci_effect)/2)) %>%
  mutate(LABEL = paste0(METRIC, ", ", SIZE, " clones"),
         AGE_DAYS = sapply(str_split(LABEL, " "), function(x) x[1])) %>%
  mutate(AGE_DAYS = factor(AGE_DAYS, levels=groups))
pal <- age_palette[groups %in% vj_alpha_bs_ci$AGE_DAYS]

# Configure series labels
vj_bs_series_labels <- vj_alpha_bs_ci %>% group_by(LABEL) %>%
  summarise(D_MAX = max(D_MEAN), D_MIN = min(D_MEAN)) %>%
  mutate(TAKE_MIN = c(TRUE, TRUE, TRUE, FALSE),
         D_SPACE = c(0.1, 0.07, 0.08, 0.06),
         D_LAB = ifelse(TAKE_MIN, D_MIN - D_SPACE,
                        D_MAX + D_SPACE))

# Make plot
g_vj_bs <- ggplot(vj_alpha_bs_ci, aes(x=Q)) +
  geom_ribbon(aes(ymin=D_LOWER_EMPIRICAL, ymax=D_UPPER_EMPIRICAL, fill=AGE_DAYS, linetype=SIZE), alpha = 0.3) +
  geom_line(aes(colour=AGE_DAYS, y=D_MEAN, linetype=SIZE),
            size = 0.3) +
  geom_text(aes(x=4, y=D_LAB, label=LABEL), 
            data = vj_bs_series_labels,
            size = fontsize_base * 5/14 * 0.6, hjust = 1, vjust = 0.5) +
  scale_colour_manual(values = pal, name = "Age group (days)") +
  scale_fill_manual(values = pal, name = "Age group (days)") +
  scale_x_continuous(name="Diversity order (q)", #expand = c(0,0),
                     limits = c(0, 4), 
                     breaks = seq(0, 4, 1)) +
  scale_y_continuous(name="α-div ratio (old/young)", expand = c(0,0),
                     limits = c(0, 1.1), 
                     breaks = seq(0, 1, 0.2)) +
  coord_fixed(ratio = 4/1.1) +
  theme_base + theme(legend.position = "none")

fig_f <- g_vj_bs

cat("done.\n")

#==============================================================================
# G - BETA SPECTRA
#==============================================================================

cat("\nGenerating VJ beta-diversity (all clones)...")

vj_beta <- import_div(vj_beta_path) %>%
  mutate(AGE_DAYS = factor(AGE_DAYS, levels = groups)) %>%
  filter(DIVTYPE == "beta") %>%
  mutate(D = (D-1)/(N_GROUP-1), D_UPPER_EMPIRICAL = (D_UPPER_EMPIRICAL-1)/(N_GROUP-1),
         D_LOWER_EMPIRICAL = (D_LOWER_EMPIRICAL-1)/(N_GROUP-1), D_SD = D_SD/(N_GROUP-1))

# Make spectrum plot
g_vj_beta <- make_spectrum_plot(vj_beta, NULL,
                                "Relative VJ beta diversity", 0.25, 0.2, 0.1)

fig_g <- g_vj_beta

cat("done.\n")

#==============================================================================
# H - RDI BOXPLOTS
#==============================================================================

cat("\nGenerating RDI boxplots...")

ctab <- import_div(ctab_path)

# Filter ambiguous segment calls
ctab_ambig <- filter(ctab, HAS_VJ, !VJ_AMBIG)

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
indivs <- ctab %>% group_by(AGE_DAYS, INDIVIDUAL) %>% summarise %>%
  ungroup() %>% arrange(INDIVIDUAL) %>% mutate(N = row_number())
indiv_index <- setNames(indivs$N, indivs$INDIVIDUAL)
indiv_ages  <-  setNames(indivs$AGE_DAYS, indivs$INDIVIDUAL)
dist_tab <- rdi %>% as.matrix %>% 
  melt(varnames = c("ID1", "ID2"), value.name = "RDI") %>% as_tibble %>%
  mutate(ROW = indiv_index[ID1], COL = indiv_index[ID2],
         AGE1 = indiv_ages[ID1], AGE2 = indiv_ages[ID2])

# Remove redundant distances and restrict to matching age groups
dist_tab_filtered <- filter(dist_tab, ROW > COL, AGE1 == AGE2) %>%
  select(AGE_DAYS = AGE1, ID1, ID2, RDI)

# Define maximum distances for annotation positioning
dist_max_tab <- dist_tab_filtered %>% group_by(AGE_DAYS) %>%
  summarise(DIST_MAX = max(RDI))

# Prepare MWU significance tables
signif_level <- 0.05
mwu_grid_all <- dist_tab_filtered %>% get_mwu_grid %>%
  filter(P <= signif_level) %>%
  inner_join(dist_max_tab %>% rename(AGE1 = AGE_DAYS), by="AGE1") %>%
  rename(DIST_MAX_AGE1 = DIST_MAX) %>%
  inner_join(dist_max_tab %>% rename(AGE2 = AGE_DAYS), by="AGE2") %>%
  rename(DIST_MAX_AGE2 = DIST_MAX) %>%
  mutate(DIST_MAX = pmax(DIST_MAX_AGE1,DIST_MAX_AGE2)) %>%
  mutate(DIST_START = min(DIST_MAX)) %>%
  select(-DIST_MAX_AGE1, -DIST_MAX_AGE2) %>%
  mutate(AGE1 = as.numeric(AGE1), AGE2 = as.numeric(AGE2),
         AGE_AVG = (AGE1 + AGE2)/2, AGE_MIN = pmin(AGE1,AGE2),
         AGE_DIFF = pmax(AGE1, AGE2) - pmin(AGE1,AGE2)) %>%
  arrange(AGE_MIN, AGE_DIFF) %>%
  mutate(Y_BAR = DIST_START + c(0.2, 0.6, 1.0),
         Y_LAB = Y_BAR + 0.1,
         LABEL = signif_stars(P))


# Make RDI boxplot
g_rdi_all <- ggplot(dist_tab_filtered, aes(x=as.numeric(AGE_DAYS), y=RDI)) +
  geom_boxplot(aes(fill = factor(AGE_DAYS, levels = groups)), 
               outlier.shape = NA) +
  geom_point(size = 0.7, alpha = 0.2, shape=16) +
  scale_fill_manual(values = age_palette, name = "Age group (days)") +
  scale_colour_manual(values = age_palette, name = "Age group (days)") +
  scale_y_continuous(name="Pairwise RDI distance", limits=c(3,6.7),
                     breaks = seq(3,6,1)) +
  scale_x_continuous(name="Age (days)", limits=c(30, 140),
                     breaks=seq(30, 120, 30)) + 
  coord_fixed((140-30)/(6.7-3)) +
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

fig_h <- g_rdi_grob

cat("done.\n")

#==============================================================================
# I - PCOA
#==============================================================================

cat("\nGenerating RDI PCoA plot...")

indivs <- ctab %>% group_by(AGE_DAYS, INDIVIDUAL) %>% summarise %>% ungroup

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
                               colour = factor(AGE_DAYS, levels = groups))) + 
  geom_point(size = 1, alpha = 0.8, shape = 16) +
  scale_colour_manual(values = age_palette, name = "Age group (days)") +
  facet_wrap(~factor(AGE_DAYS, levels = groups), scales = "fixed") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=axis_width) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=axis_width) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(name=paste0("Principle co-ord. 1 (", pc_var[1], "%)"),
                     limits = c(-pcoa_scale,pcoa_scale), breaks = seq(-2,2,2)) +
  scale_y_continuous(name=paste0("Principle co-ord. 2 (", pc_var[2], "%)"),
                     limits = c(-pcoa_scale,pcoa_scale), breaks = seq(-2,2,2)) +
  theme_base +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

fig_i <- g_pcoa

cat("done.\n")

#==============================================================================
# LEGEND
#==============================================================================

cat("\nGenerating legend...")

g_legend <- clone_alpha %>% filter(Q %in% c(0,1,2)) %>%
  ggplot(aes(x=Q)) +
  geom_ribbon(aes(ymin=D_LOWER_EMPIRICAL, ymax=D_UPPER_EMPIRICAL, fill=AGE_DAYS),
              alpha = 0.4) +
  geom_line(aes(colour=AGE_DAYS, y=D), size = 0.5) +
  geom_point(aes(colour=AGE_DAYS, y=D), size = 1.3) +
  scale_colour_manual(values = age_palette, name = "Age group (days)") +
  scale_fill_manual(values = age_palette, name = "Age group (days)") +
  guides(colour = guide_legend(keyheight=legend_keysize, keywidth = legend_keysize),
         fill   = guide_legend(keyheight=legend_keysize, keywidth = legend_keysize)) +
  theme_base + theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(margin = margin(r=4))
    )
legend_out <- get_legend(g_legend)

#==============================================================================
# MAKE FIGURE GRID
#==============================================================================

cat("\nPreparing output figure...")

fig_row1 <- grid_fig(fig_a, ncol = 1, labels="A", vjust = 1)
fig_row2 <- grid_fig(fig_b, fig_c, ncol = 2, labels = c("B", "C"))
fig_row3 <- grid_fig(fig_d, fig_e, fig_f, ncol = 3, labels = c("D", "E", "F"))
fig_row4 <- grid_fig(fig_g, fig_h, fig_i, ncol = 3, labels = c("G", "H", "I"))
fig <- grid_fig(fig_row1, fig_row2, fig_row3, fig_row4, nrow = 4, 
                 labels = NULL, rel_heights = c(0.75,1,1,1))
fig_legend <- plot_grid(fig, legend_out, nrow = 2, rel_heights = c(1, 0.03))

cat("done.\n")

#==============================================================================
# SAVE OUTPUT
#==============================================================================

cat("\nSaving output...")

save_fig(snakemake@output[[1]], fig_legend, plot_height, plot_width)

cat("done.\n")