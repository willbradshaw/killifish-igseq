###############################################################################
## SUPPLEMENTARY FIGURE                                                      ##
## Clone-size distributions from whole-body ageing samples                  ##
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
ctab_path    <- snakemake@input[[1]]

# Define output parameters
out_path <- snakemake@output[[1]]
out_width <- snakemake@params[["plot_width"]]
out_height <- snakemake@params[["plot_ratio"]] * out_width

# Modify theme
theme_base <- theme_base + theme(
  legend.text = element_text(margin = margin(r = 0.2, unit = "cm"))
)

# Configure run parameters
palette <- ageing_palette
age_groups <- ageing_groups
size_limit <- 5 # Size threshold for large clones

cat("done.\n")

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

cat("\nImporting data...")
ctab <- suppressMessages(read_tsv(ctab_path))
cat("done.\n")

#------------------------------------------------------------------------------
# EVALUATE CLONE SIZES
#------------------------------------------------------------------------------

cat("\nProcessing data...")

# Count unique sequences per clone
clones <- ctab %>% filter(!is.na(CLONE)) %>%
  group_by(AGE_DAYS, INDIVIDUAL, CLONE, SEQUENCE_INPUT, BEST_VJ_CALL) %>%
  summarise() %>% group_by(AGE_DAYS, INDIVIDUAL, CLONE, BEST_VJ_CALL) %>%
  summarise(NSEQ = n()) %>% mutate(LARGE = NSEQ >= 5)

# Get clone size distribution for each individual
clone_sizes <- clones %>% group_by(AGE_DAYS, INDIVIDUAL, CLONE) %>%
  summarise(NSEQ = sum(NSEQ)) %>% group_by(AGE_DAYS, INDIVIDUAL, NSEQ) %>%
  summarise(NCLONE = n()) %>%
  mutate(NCLONE_PC = NCLONE/sum(NCLONE)*100, NCLONE_PC_CM = cumsum(NCLONE_PC))

# Get clone sequence distribution for each individual
clone_seqs <- clone_sizes %>%
  mutate(NSEQ_TOTAL = NSEQ * NCLONE, 
         NSEQ_TOTAL_PC = NSEQ_TOTAL/sum(NSEQ_TOTAL) * 100,
         NSEQ_TOTAL_PC_CM = cumsum(NSEQ_TOTAL_PC))

# Get V/J size distributions
clones_vj <- clones %>% group_by(AGE_DAYS, INDIVIDUAL, BEST_VJ_CALL, LARGE) %>%
  summarise(NSEQ = sum(NSEQ), NCLONE = n())

clones_avg <- clones_vj %>% group_by(AGE_DAYS, INDIVIDUAL, LARGE) %>%
  summarise(NSEQ = sum(NSEQ), NCLONE = sum(NCLONE)) %>%
  mutate(NCLONE_TOTAL = sum(NCLONE), NSEQ_TOTAL = sum(NSEQ),
         NSEQ_PC = NSEQ/sum(NSEQ), NCLONE_PC = NCLONE/sum(NCLONE)) %>%
  ungroup() %>% filter(LARGE) %>% select(-LARGE)

seqs_vj_prop <- clones_vj %>% group_by(AGE_DAYS, INDIVIDUAL, BEST_VJ_CALL) %>%
  summarise(NSEQ = sum(NSEQ), NCLONE = sum(NCLONE)) %>%
  mutate(SEQS_PER_CLONE = NSEQ/NCLONE)

seqs_vj_group <- seqs_vj_prop %>%
  mutate(NSEQ_OOM = floor(log10(NSEQ)), LARGE = SEQS_PER_CLONE > size_limit) %>%
  group_by(AGE_DAYS, NSEQ_OOM, LARGE) %>% summarise(N = n()) %>%
  group_by(AGE_DAYS, NSEQ_OOM) %>% mutate(PC = N/sum(N) * 100) %>%
  mutate(X = 10^(NSEQ_OOM + 0.05), Y = ifelse(LARGE, 100, size_limit * 0.95))

cat("done.\n")

#------------------------------------------------------------------------------
# PLOT CLONE SIZE STATISTICS
#------------------------------------------------------------------------------

cat("\nMaking plots...")

# Clone size distributions
g_clone_sizes <- ggplot(clone_sizes) +
  geom_line(aes(x=NSEQ, y=NCLONE_PC_CM, group=INDIVIDUAL,
                colour = factor(AGE_DAYS, levels = age_groups))) +
  geom_vline(xintercept = 5, linetype = "dashed") +
  #facet_wrap(~AGE_DAYS, scales="free") +
  scale_x_log10(name = "Unique sequences per clone",
                breaks = c(1,5,10,100,1000), limits=c(1, 1000),
                expand = c(0,0)) +
  scale_y_continuous(name="Cumulative % of clones",
                     limits=c(60,101), expand = c(0,0)) +
  coord_fixed(ratio = (log10(1000)-log10(1))/(101-60)) +
  scale_colour_manual(values = palette, name = "Age (days)") +
  theme_base + theme(legend.position = "none")


# V/J distribution
g_v_seqs <- ggplot(seqs_vj_prop) +
  geom_point(aes(x=NSEQ, y=SEQS_PER_CLONE, 
                 colour = factor(AGE_DAYS, levels = age_groups)), 
             alpha = 0.5, size = 1, shape = 16) + 
  #facet_wrap(~AGE_DAYS, scales="fixed") +
  geom_hline(yintercept = size_limit, linetype = "dashed") +
  scale_colour_manual(values = palette, name = "Age (days)") +
  scale_x_log10(name = "Unique sequences per V/J", limits=c(1,1000),
                expand = c(0,0)) +
  scale_y_log10(name = "Mean sequences per clone", limits=c(1,200),
                breaks = c(1,size_limit,10,100), expand=c(0,0)) +
  coord_fixed(ratio = (log10(1000))/(log10(200))) +
  theme_base + theme(legend.position = "none")

#------------------------------------------------------------------------------
# MAKE COMBO LEGEND
#------------------------------------------------------------------------------

keysize <- 0.7

g_legend <- ggplot(clone_sizes) +
  geom_line(aes(x=NSEQ, y=NCLONE_PC_CM, group=INDIVIDUAL,
                colour = factor(AGE_DAYS, levels = age_groups))) +
  geom_point(aes(x=NSEQ, y=NCLONE_PC_CM, group=INDIVIDUAL,
                colour = factor(AGE_DAYS, levels = age_groups)),
             size = 1, shape = 16) +
  scale_colour_manual(values = palette, name = "Age (days)") +
  guides(colour = guide_legend(keyheight=keysize, keywidth = keysize),
         fill   = guide_legend(keyheight=keysize, keywidth = keysize)) +
  theme_base + theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(margin = margin(r=4))
  )
legend_out <- get_legend(g_legend)

#------------------------------------------------------------------------------
# MAKE GRID PLOT
#------------------------------------------------------------------------------

fig <- grid_fig(g_clone_sizes, g_v_seqs, ncol = 2, vjust = 1)
fig_out <- plot_grid(fig, legend_out, nrow = 2, rel_heights = c(1, 0.08))

cat("done.\n")

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

cat("\nSaving output...")
save_fig(out_path, fig_out, out_height, out_width)
cat("...done.\n")