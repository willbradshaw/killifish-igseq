###############################################################################
## FIGURE                                                                    ##
## GLM FIGURE MAIN                                                           ##
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

# Set parameters
palette_exp <- c("#1F78B4", "#E31A1C")
qvals_short <- c(0,1,2)
qvals_long <- c(0,1,1.5,2,3,4)
individuals_excluded <- snakemake@params[["individuals_excluded"]]

# Configure input paths
rarefy_path <- snakemake@input[["rarefy"]]
gut_clone_path <- snakemake@input[["gut_clone"]]
gut_vj_all_path <- snakemake@input[["gut_vj_all"]]
gut_vj_large_path <- snakemake@input[["gut_vj_large"]]
gut_vj_small_path <- snakemake@input[["gut_vj_small"]]
age_clone_path <- snakemake@input[["age_clone"]]
age_vj_all_path <- snakemake@input[["age_vj_all"]]
age_vj_large_path <- snakemake@input[["age_vj_large"]]
age_vj_small_path <- snakemake@input[["age_vj_small"]]

#==============================================================================
# AUXILIARY FUNCTIONS
#==============================================================================

get_glm <- function(tab, q, type, family = Gamma(),
                    fl = formula("D_NORM~AGE_DAYS*SAMPLE")){
  tab %>% filter(Q == q, TYPE == type) %>% 
    {if (is.na(family)[1]) lm(fl, data = .)
      else glm(fl, family = family, data = .)} 
}

make_glm_plot <- function(tab = tab_norm, predictions = glm_predictions,
                          pvals = glm_sample_age_pvalues,
                          qvals = qvals, types = spectrum_types){
  tab %>% ungroup %>% filter(Q %in% qvals, TYPE %in% types) %>% 
    mutate(TYPE = factor(TYPE, levels = types)) %>%
    arrange(TYPE) %>%
    ggplot(aes(x=AGE_DAYS)) +
    geom_line(aes(y=D_NORM,colour=SAMPLE), 
              data=predictions %>% filter(Q %in% qvals, TYPE %in% types)) +
    geom_boxplot(aes(y=D_NORM,fill=SAMPLE, group=interaction(SAMPLE,AGE_DAYS)),
                 outlier.shape = NA) +
    geom_point(aes(y=D_NORM), size = 0.2, alpha = 0.5, shape=16) +
    geom_linerange(aes(x=X_LINE,ymin=Y_MIN,ymax=Y_MAX), size = 0.2,
                   data=pvals %>% filter(Q %in% qvals, TYPE %in% types))+
    geom_text(aes(x=X_LABEL, y=Y_LABEL, label=LABEL),
              angle=270, size=2, vjust=1, hjust=0.5, family=font,
              data = pvals %>% filter(Q %in% qvals, TYPE %in% types)) +
    facet_grid(factor(TYPE, levels = types) ~ Q, scales = "free",
               labeller = labeller(Q = function(x) paste("Diversity order", x))) +
    scale_y_continuous(name = "Relative diversity", breaks = seq(0,2),
                       limits = c(0,2.1)) +
    scale_x_continuous(name = "Age (days)", breaks = seq(30,140,30),
                       limits = c(30,150)) +
    scale_fill_manual(name = "Sample type", values=palette_exp) +
    scale_colour_manual(name = "Sample type", values=palette_exp) +
    theme_base +
    theme(plot.margin = margin(0,0,0,5,"mm"),
          panel.spacing = unit(2, "mm"),
          strip.text.y = element_text(face = "bold"),
          strip.text.x = element_text(face = "bold"))
}

cat("...done.\n")

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

cat("\nImporting data...")

# Rarefied clones
rarefy <- import_div(rarefy_path) %>%
  mutate(SAMPLE = ifelse(EXPERIMENT == "gut", "Gut", "Whole body"))

# Gut spectra
gut_clone <- import_div(gut_clone_path) %>% 
  mutate(SAMPLE = "Gut", TYPE = "Clonal repertoire")
gut_vj_all <- import_div(gut_vj_all_path) %>% 
  mutate(SAMPLE = "Gut", TYPE = "V/J repertoire")
gut_vj_large <- import_div(gut_vj_large_path) %>% 
  mutate(SAMPLE = "Gut", TYPE = "V/J (large clones)")
gut_vj_small <- import_div(gut_vj_small_path) %>% 
  mutate(SAMPLE = "Gut", TYPE = "V/J (small clones)")
div_gut <- bind_rows(gut_clone, gut_vj_all, gut_vj_large, gut_vj_small)

# Ageing
age_clone <- import_div(age_clone_path) %>% 
  mutate(SAMPLE = "Whole body", TYPE = "Clonal repertoire")
age_vj_all <- import_div(age_vj_all_path) %>% 
  mutate(SAMPLE = "Whole body", TYPE = "V/J repertoire")
age_vj_large <- import_div(age_vj_large_path) %>% 
  mutate(SAMPLE = "Whole body", TYPE = "V/J (large clones)")
age_vj_small <- import_div(age_vj_small_path) %>% 
  mutate(SAMPLE = "Whole body", TYPE = "V/J (small clones)")
div_age <- bind_rows(age_clone, age_vj_all, age_vj_large, age_vj_small)

cat("...done.\n")

#------------------------------------------------------------------------------
# SUBFIGURE A (RAREFACTION)
#------------------------------------------------------------------------------

cat("\nGenerating rarefaction plots...")

# Small clones
fig_a_small <- rarefy %>% filter(METRIC == "N_CLONES_SMALL", 
                            !INDIVIDUAL %in% individuals_excluded) %>%
  ggplot +
  geom_line(aes(x=SAMPLE_SIZE, y=MEAN, colour = SAMPLE,
                group = interaction(SAMPLE,INDIVIDUAL))) +
  scale_color_manual(values = palette_exp, name = "Sample type") +
  scale_x_continuous(breaks = seq(0, 10000, 2000),
                     labels = function(x) x/1000,
                     name = "Number of UMI groups (×1000)") +
  scale_y_continuous(labels = function(y) y/1000,
                     name = "No. of small clones (×1000)") +
  geom_ribbon(aes(x=SAMPLE_SIZE, ymin=MEAN-SD, ymax = MEAN + SD,
                  fill = SAMPLE,
                  group = interaction(SAMPLE,INDIVIDUAL)), alpha = 0.4) +
  scale_fill_manual(values = palette_exp, name = "Sample type") +
  theme_classic() + theme_base +
  theme(legend.key.size = unit(0.3, "cm"),
        aspect.ratio = 1.1,
        plot.margin = margin(b = 0, l = 0.1, r = 0.1, unit = "cm"))

# Large clones
fig_a_large <- rarefy %>% filter(METRIC == "N_CLONES_LARGE", 
                            !INDIVIDUAL %in% individuals_excluded) %>%
  ggplot +
  geom_line(aes(x=SAMPLE_SIZE, y=MEAN, colour = SAMPLE,
                group = interaction(SAMPLE,INDIVIDUAL))) +
  scale_color_manual(values = palette_exp, name = "Sample type") +
  scale_x_continuous(breaks = seq(0, 10000, 2000),
                     labels = function(x) x/1000,
                     name = "Number of UMI groups (×1000)") +
  scale_y_continuous(labels = function(y) y/100,
                     name = "No. of large clones (×100)") +
  geom_ribbon(aes(x=SAMPLE_SIZE, ymin=MEAN-SD, ymax = MEAN + SD,
                  fill = SAMPLE,
                  group = interaction(SAMPLE,INDIVIDUAL)), alpha = 0.4) +
  scale_fill_manual(values = palette_exp, name = "Sample type") +
  theme_classic() + theme_base +
  theme(legend.key.size = unit(0.3, "cm"),
        aspect.ratio = 1.1,
        plot.margin = margin(b = 0, l = 0.1, r = 0.1, unit = "cm"))

# Percent large clones
fig_a_pc <- rarefy %>% filter(METRIC == "PC_CLONES_SMALL", 
                            !INDIVIDUAL %in% individuals_excluded) %>%
  ggplot +
  geom_line(aes(x=SAMPLE_SIZE, y=1-MEAN, colour = SAMPLE,
                group = interaction(SAMPLE,INDIVIDUAL))) +
  scale_color_manual(values = palette_exp, name = "Sample type") +
  scale_x_continuous(breaks = seq(0, 10000, 2000),
                     labels = function(x) x/1000,
                     name = "Number of UMI groups (×1000)") +
  scale_y_continuous(labels = function(y) y*100, breaks = seq(0, 0.6, 0.2),
                     name = "% of large clones", limits = c(0, 0.7)) +
  geom_ribbon(aes(x=SAMPLE_SIZE, ymin=1-MEAN-SD, ymax = 1-MEAN + SD,
                  fill = SAMPLE,
                  group = interaction(SAMPLE,INDIVIDUAL)), alpha = 0.4) +
  scale_fill_manual(values = palette_exp, name = "Sample type") +
  theme_classic() + theme_base +
  theme(legend.key.size = unit(0.3, "cm"),
        aspect.ratio = 1.1,
        plot.margin = margin(b = 0, l = 0.1, r = 0.1, unit = "cm"))

rarefy_legend <- get_legend(fig_a_pc + theme(
  legend.justification = "center"
))

fig_a_raw <- plot_grid(fig_a_small + theme(
  legend.position = "none", axis.title.x = element_blank()),
  fig_a_large + theme(legend.position = "none"),
  fig_a_pc + theme(legend.position = "none", 
                   axis.title.x = element_blank()),
  nrow = 1, labels = NULL, align = "hv", axis="l")

fig_a <- plot_grid(fig_a_raw, rarefy_legend, nrow = 2,
                   rel_heights = c(1, 0.06), labels = NULL)

cat("done.\n")

#------------------------------------------------------------------------------
# SUBFIGURE B (GLMs)
#------------------------------------------------------------------------------

cat("\nGenerating GLM plots...")

# Combine spectra together
tab_b <- bind_rows(div_gut, div_age) %>%
  mutate(AGE_DAYS = as.numeric(AGE_DAYS), AGE_WEEKS = as.numeric(AGE_WEEKS),
         AGE_DAYS = ifelse(is.na(AGE_DAYS), AGE_WEEKS*7, AGE_DAYS),
         AGE_WEEKS = ifelse(is.na(AGE_WEEKS), AGE_DAYS/7, AGE_WEEKS)) %>%
  select(-matches("^D_"), -matches("^E"), -DIVTYPE, -N_GROUP)

# Normalise diversities by youngest group
tab_norm <- tab_b %>% group_by(TYPE, SAMPLE, Q) %>%
  filter(AGE_DAYS == min(AGE_DAYS)) %>% summarise(D_REF = mean(D)) %>%
  full_join(tab_b, by=c("TYPE", "SAMPLE", "Q")) %>% 
  mutate(D_NORM = D/D_REF)

# Get GLMs for each spectrum type and Q-value
spectrum_types <- c("Clonal repertoire", "V/J repertoire", "V/J (large clones)",
                    "V/J (small clones)")
glm_list <- lapply(spectrum_types, function(type)
  lapply(qvals_long, function(q) get_glm(tab_norm, q, type)))

# Get GLM predictions for each model
ages_gut <- tab_b %>% filter(SAMPLE == "Gut") %>% pull("AGE_DAYS") %>% 
  range %>% (function(x) tibble(AGE_DAYS = seq(x[1],x[2]), SAMPLE = "Gut"))
ages_age <- tab_b %>% filter(SAMPLE == "Whole body") %>% pull("AGE_DAYS") %>% 
  range %>% (function(x) tibble(AGE_DAYS = seq(x[1],x[2]), SAMPLE = "Whole body"))
glm_predictions_gut <- lapply(glm_list, function(x) lapply(x, function(g)
  ages_gut %>% mutate(PREDICTION = predict.glm(g, ages_gut, type = "response")))) %>%
  melt(id.var=c("SAMPLE","AGE_DAYS"), value.name = "PREDICTION") %>% as_tibble
glm_predictions_age <- lapply(glm_list, function(x) lapply(x, function(g)
  ages_age %>% mutate(PREDICTION = predict.glm(g, ages_age, type = "response")))) %>%
  melt(id.var=c("SAMPLE","AGE_DAYS"), value.name = "PREDICTION") %>% as_tibble
glm_predictions <- bind_rows(glm_predictions_gut, glm_predictions_age) %>%
  mutate(TYPE = spectrum_types[L1], Q = qvals_long[L2]) %>% 
  rename(D_NORM = PREDICTION) %>% select(-L1,-L2,-variable)

# Test for a significant sample-type effect on age-sensitivity
glm_sample_age_pvalues <- lapply(glm_list, function(x) lapply(x, function(g)
  summary(g)$coefficients[4,4])) %>% melt(value.name = "P") %>% as_tibble %>%
  mutate(TYPE = spectrum_types[L1], Q = qvals_long[L2]) %>% select(-L1,-L2) %>%
  mutate(LABEL = signif_stars(P))

# Get positions of significance labels for plots
label_coordinates <- glm_predictions %>% group_by(SAMPLE) %>%
  filter(AGE_DAYS == max(AGE_DAYS)) %>% 
  spread(SAMPLE,D_NORM,fill=-Inf,sep="_") %>%
  rename(Y_GUT = SAMPLE_Gut, Y_AGE = `SAMPLE_Whole body`) %>%
  group_by(TYPE, Q) %>% 
  summarise(X = max(AGE_DAYS), Y_GUT = max(Y_GUT), Y_AGE = max(Y_AGE)) %>%
  mutate(Y_MIN = pmin(Y_GUT, Y_AGE), Y_MAX = pmax(Y_GUT, Y_AGE),
         Y_LABEL = (Y_GUT + Y_AGE)/2,
         X_LINE = round(X/10)*10 + 5, X_LABEL = X_LINE + 10)

glm_sample_age_pvalues <- full_join(glm_sample_age_pvalues, label_coordinates,
                                    by=c("TYPE", "Q")) %>%
  mutate(LABEL = signif_stars(P))

# Make plots
fig_b <- make_glm_plot(qvals = qvals_short, types = spectrum_types[-4]) +
  theme(legend.key.size = unit(0.3, "cm"),
        plot.margin = margin(t=0.5, r = 0.1, l = 0.1, b = 0, unit = "cm"),
        legend.margin = margin(t=0, b=0.1))
fig_si <- make_glm_plot(qvals = qvals_long, types = spectrum_types) +
  theme(legend.key.size = unit(0.3, "cm"))

cat("done.\n")

#------------------------------------------------------------------------------
# PUT IT ALL TOGETHER
#------------------------------------------------------------------------------

cat("\nAssembling output figure...")
fig_all <- plot_grid(fig_a, fig_b, ncol = 1, nrow = 2, labels = LETTERS[1:2],
                     label_size = fontsize_base * fontscale_label,
                     label_fontfamily = titlefont, label_colour = "black",
                     label_fontface = "plain",
                     rel_heights = c(1,2))
cat("done.\n")

cat("\nSaving output...")
save_fig(snakemake@output[["main"]], fig_all,
         snakemake@params[["plot_width_main"]] * snakemake@params[["plot_ratio_main"]],
         snakemake@params[["plot_width_main"]])

save_fig(snakemake@output[["si"]], fig_si,
         snakemake@params[["plot_width_si"]] * snakemake@params[["plot_ratio_si"]],
         snakemake@params[["plot_width_si"]])
cat("done.\n")
