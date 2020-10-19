###############################################################################
## MAIN TEXT FIGURE                                                          ##
## Repertoire composition and diversity in whole-body killifish samples      ##
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
plot_width <- snakemake@params[["plot_width"]]
plot_height <- plot_width * snakemake@params[["plot_ratio"]] 

# Specify input paths
entropies_path <- snakemake@input[["entropies_pilot_group"]]
indels_path    <- snakemake@input[["indels"]]
ageing_path    <- snakemake@input[["entropies_ageing_solo"]]

cat("done.\n")

#==============================================================================
# IMPORT DATA
#==============================================================================

cat("\nImporting data...")

# Pilot group entropies
entropies <- import_entropies(entropies_path)

# Pilot group indels
indels <- import_div(indels_path,
                     cols(p = "d", n = "i", .default = "c"))

# Ageing individual entropies
ageing <- import_entropies(ageing_path) %>%
    mutate(age_group = as.integer(sub("-.*", "", id)),
           age_days = as.integer(ageing_groups[age_group]))

cat("done.\n")

#==============================================================================
# A - ENTROPY GRID
#==============================================================================

cat("\nGenerating entropy grid...")

# Process entropies for plotting
tab_a <- entropy_table(entropies)

# Plot entropy composition
fig_a <- plot_entropy_table(tab_a) + theme(
    plot.margin = margin(b=0.1, t=0.1, unit="cm"),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.background = element_blank(),panel.background = element_blank(),
    axis.title.x = element_blank(), axis.title.x.top = element_blank(),
    axis.title.x.bottom = element_blank(), axis.title.y = element_blank(),
    legend.margin = margin(0), strip.text = element_blank()
    )

cat("done.\n")

#==============================================================================
# B&C - INDELS
#==============================================================================

cat("\nGenerating indel plots...")

theme_indel = theme_base + theme(
    plot.margin = margin(0.1,0.1,0,0.2, unit="cm"),
    legend.margin = margin(0,0,0.2,0, unit="cm"),
    legend.title = element_text(margin=margin(r=0.1, unit="cm")),
    axis.title.x = element_text(margin = margin(2,0,0,0)),
    axis.title.y = element_text(margin = margin(0,2,0,0)),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    aspect.ratio = 1,
    )

fig_b <- ggplot(mapping = aes(x=n, y=p, colour=event)) +
    geom_line(data = indels %>% filter(event == "vd_ins")) +
    geom_line(data = indels %>% filter(event == "dj_ins")) +
    
    scale_colour_manual(values = palette_indels, name = NULL,
                        labels = c("D/J insertions", "V/D insertions")) +
    scale_y_continuous(name = "Probability (%)",
                       labels = function(y) y * 100, expand = c(0,0)) +
    scale_x_continuous(name = "# Insertions", limits = c(NA,15)) +
    theme_classic() + theme_indel + theme(
        legend.justification = c(1,1),
        legend.background = element_rect(fill=alpha("white", 0.5)),
        legend.text = element_text(margin = margin(t=0.8, b=0.8, unit="mm")),
        legend.position = c(1,1),
        # legend.spacing.y = unit(vspace, "mm"),
    )

fig_c <- ggplot(mapping = aes(x=n, y=p, colour=event)) +
    geom_line(data = indels %>% filter(event == "v_3_del")) +
    geom_line(data = indels %>% filter(event == "j_5_del")) +
    geom_line(data = indels %>% filter(event == "d_3_del")) +
    geom_line(data = indels %>% filter(event == "d_5_del")) +
    scale_colour_manual(values = palette_indels, name = NULL,
                        labels = c("5' D-deletions", "3' D-deletions",
                                   "5' J-deletions", "3' V-deletions")
                        ) +
    scale_y_continuous(name = "Probability (%)",
                       labels = function(y) y * 100, expand = c(0,0),
                       limits = c(0,0.25)) +
    scale_x_continuous(name = "# Deletions", limits = c(NA,27)) +
    theme_classic() + theme_indel + theme(
        legend.justification = c(1,1),
        legend.background = element_rect(fill=alpha("white", 0.5)),
        legend.text = element_text(margin = margin(t=0.7, b=0.7, unit="mm")),
        legend.position = c(1,1),
        # legend.spacing.y = unit(vspace, "mm"),
    )

cat("done.\n")

#==============================================================================
# D - AGE EFFECT ON ENTROPY
#==============================================================================

cat("\nGenerating age/entropy plot...")

# Separate entropies by category
tab_d <- ageing %>% filter(short == "total") %>% mutate(d = 2^h)

# Test for an age effect
p_d <- kruskal.test(formula = d ~ age_days, data=tab_d)$p.value

# Make plot
fig_d_raw <- ggplot(tab_d) +
    geom_boxplot(aes(x=age_days, y=h, fill=factor(age_days, levels=ageing_groups),
                     group=age_days), outlier.shape = NA) +
    geom_point(aes(x=age_days, y=h), size = 0.5, alpha = 0.5, shape = 16) +
    geom_text(x=128, y=23, label=paste0("p = ", round(p_d, 2)), family = font,
              size = fontsize_base * 5/14, hjust=1, vjust=0, face = "plain") +
    scale_fill_manual(values = ageing_palette, name="Age (days)") +
    scale_x_continuous(name = "Age (days)") +
    scale_y_continuous(name = "Generative diversity (bits)") +
    theme_classic() + theme_base + theme(
         plot.margin = margin(0.1,0.1,0,0.2, unit="cm"),
         legend.position = "none",
         axis.title.x = element_text(margin = margin(2,0,0,0)),
         axis.title.y = element_text(margin = margin(0,2,0,0)),
         panel.background = element_blank(),
         plot.background = element_blank(),
         aspect.ratio = 1,
         )

# Remove boxplot borders
fig_d_grob <- ggplotGrob(fig_d_raw)
for (n in 1:length(fig_d_grob$grobs[[6]]$children[[3]]$children)){
    h <- which(grepl("geom_crossbar",names(fig_d_grob$grobs[[6]]$children[[3]]$children[[n]]$children)))
    fig_d_grob$grobs[[6]]$children[[3]]$children[[n]]$children[[h]]$children[[1]]$gp$col <- fig_d_grob$grobs[[6]]$children[[3]]$children[[n]]$children[[h]]$children[[1]]$gp$fill
}

fig_d <- fig_d_grob

cat("done.\n")

#==============================================================================
# MAKE FIGURE GRID
#==============================================================================

cat("\nPreparing output figure...")

fig_row1 <- grid_fig(fig_a, ncol = 1, labels="A")
fig_row2 <- grid_fig(fig_b, fig_c, fig_d, ncol = 3, labels = c("B", "C", "D"))
fig <- grid_fig(fig_row1, fig_row2, nrow = 2,
                 labels = NULL, rel_heights = c(0.55,1))

cat("done.\n")

#==============================================================================
# SAVE OUTPUT
#==============================================================================

cat("\nSaving PNG output...")
save_fig(snakemake@output[["png"]], fig, plot_height, plot_width, device = "png")
cat("done.\n")
cat("\nSaving RDS output...")
saveRDS(fig, snakemake@output[["rds"]])
cat("done.\n")
