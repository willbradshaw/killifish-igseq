###############################################################################
## SUPPLEMENTARY FIGURE                                                      ##
## Cross-replicate RDI plots for pilot dataset                               ##
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
ctab_path    <- snakemake@input[["changeo_db"]]
repo_path <- snakemake@input[["repo_dir"]]
lib_path  <- snakemake@input[["lib_dir"]]
# ctab_path <- "../data_processed/changeo/pilot-seqs-all.tab.gz"
# repo_path <- "env/minicran_rdi_repo/"
# lib_path  <- "env/minicran_rdi_lib/"

# Install RDI from minicran
repo <- normalizePath(repo_path) %>% paste0("file://", .)
lib <- lib_path
suppressMessages(suppressWarnings(install.packages("rdi", repos = repo,
                                                   lib = lib)))
library(rdi, lib.loc = lib)


# Define output parameters
out_path <- snakemake@output[[1]]
out_width <- snakemake@params[["plot_width"]]
out_height <- snakemake@params[["plot_ratio"]] * out_width
# out_path <- "output_files/dev_si_pilot-rdi.png"
# out_width <- 11
# out_ratio <- 0.6
# out_height <- out_width * out_ratio

# Define run parameters
rdi_distance <- "euclidean"
rdi_iterations <- 100
rdi_constant_scale <- TRUE
rdi_transform <- TRUE
clust_method <- "average"
rootedge_length <- -0.2
treeline_width <- 0.9
tiplab_offset <- -0.1
col_grey <- "#444444"
palette_tree <- c(pilot_palette, col_grey)

# Modify theme
theme_base <- theme_base + theme(
  legend.text = element_text(margin = margin(r = 0.1, unit = "cm"))
)

#==============================================================================
# Auxiliary functions
#==============================================================================


# Check tip labelling and find ancestry patterns
grepl_tips <- function(pattern, tree_object){
  # Search for a string pattern in tip labels and return all matching tips
  tree_object$tip.label[grepl(pattern, tree_object$tip.label)]
}

mrca_vec <- function(tree_object, tips){
  # Find the MRCA of a vector of tip IDs
  if (length(tips) <= 1) return(NA)
  # Convert to numeric IDs
  if (!is.integer(tips)) tips <- nodeid(tree_object, tips)
  # Iterate through nodes to find overall MRCA
  node1 <- tips[1]
  node2 <- tips[2]
  mrca <- MRCA(tree_object, node1, node2)
  tips <- tips[! tips %in% offspring(tree_object, mrca)]
  while (length(tips) > 0){
    mrca <- MRCA(tree_object, mrca, tips[1])
    tips <- tips[! tips %in% offspring(tree_object, mrca)]
  }
  return(mrca)
}

get_mrca <- function(pattern, tree_object){
  # Find the MRCA of all tips matching a string pattern
  nodes <- grepl_tips(pattern, tree_object)
  return(mrca_vec(tree_object, nodes))
}

cat("done.\n")

#==============================================================================
# Import data
#==============================================================================

cat("\nImporting data...")
ctab <- suppressMessages(read_tsv(ctab_path))
cat("done.\n")

#==============================================================================
# Process data
#==============================================================================

cat("\nPreparing RDI measurements...")

# Filter ambiguous segment calls
ctab_filtered <- filter(ctab, HAS_VJ, !VJ_AMBIG)

# Rename replicates
ctab_filtered <- mutate(ctab_filtered,
                        REPLICATE = sub("bio", "-1", REPLICATE),
                        REPLICATE = sub("lib", "-2", REPLICATE),
                        REPLICATE = sub("orig", "-3", REPLICATE))

# Compute RDI segment calls
genes <- pull(ctab_filtered, BEST_VJ_CALL)
annots <- as.character(pull(ctab_filtered, REPLICATE))
counts <- calcVDJcounts(genes = genes, seqAnnot = annots,
                        simplifyNames = TRUE, splitCommas = FALSE)

# Compute RDI distance matrix
rdi <- calcRDI(counts, subsample = TRUE,
               distMethod = rdi_distance, nIter = rdi_iterations,
               constScale = rdi_constant_scale,
               units = ifelse(rdi_transform, "lfc", "pct"))

# Perform hierarchical clustering
clust <- hclust(rdi, method = clust_method)

cat("done.\n")

#==============================================================================
# Generate dendrogram 
#==============================================================================

cat("\nGenerating dendrogram...")

# Convert to phylotibble and extract individual/replicate info
phylo <- as.phylo(clust)
phylo_tbl <- as_tibble(phylo) %>%
  mutate(individual = sub("(2-0\\d)(.*)", "\\1", label),
         replicate = sub("(2-0\\d)(.*)", "\\2", label),
         label = sub("(2-0\\d)(.*)", "\\1 \\2", label))
class(phylo_tbl) <- c("tbl_tree", class(phylo_tbl))

# Propagate individual info up tree
for (i in unique(phylo_tbl %>% filter(!is.na(individual)) %>% 
                 pull(individual))){
  ca <- get_mrca(i, phylo)
  of <- offspring(phylo_tbl, ca)
  phylo_tbl[phylo_tbl$node %in% c(ca, of$node),]$individual <- i
}
phylo_tbl[is.na(phylo_tbl$individual),]$individual <- "GR"

# Convert back to treedata
treedata_rdi <- as.treedata(phylo_tbl)
treedata_rdi@phylo$root.edge <- rootedge_length

# Make basic plot and colour by isotype
theme_tre <- theme_minimal() + theme_base +
  theme(legend.position = "none",
        legend.margin = margin(b=0,t=0.5, unit="cm"),
        plot.margin = margin(t=0, unit="cm"),
        axis.title.x = element_text(margin = margin(t=0, b=0.5, unit="cm")),
        axis.title.y = element_text(margin = margin(l = 2, t = 2, unit = "cm")))


treeplot <- revts(ggtree(treedata_rdi, aes(colour = individual), 
                         size = treeline_width)) +
  geom_rootedge(colour = col_grey, size = treeline_width) +
  scale_colour_manual(values = palette_tree) + 
  coord_flip() + 
  scale_x_reverse(labels = function(x) -x, name = "VJ-RDI", limits = c(1.5, NA),
                  breaks = seq(0,-5)) +
  scale_y_continuous(breaks=NULL, name = "Replicate") +
  geom_tiplab(size = fontsize_base * 5/14, offset = tiplab_offset, 
              family = font, angle = 270) +
  theme_minimal() + theme_base + theme(
    legend.position = "none",
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.y = element_line("#EEEEEE"),
    axis.title.x = element_text(margin = margin(t = -1.3, unit = "cm"))
    )

cat("done.\n")

#==============================================================================
# Generate PCoA plot 
#==============================================================================

cat("\nGenerating PCoA plot...")

# Perform PCoA
pc <- pcoa(rdi)

# Extract co-ordinates in first two axes into tibble for plotting
pc_tab <- tibble(REPLICATE = rownames(as.matrix(rdi)),
                 PCO1 = pc$vectors[,1],
                 PCO2 = pc$vectors[,2]) %>%
  mutate(INDIVIDUAL = sub("(2-0\\d)(.*)", "\\1", REPLICATE),
         REP = sub("(2-0\\d)(.*)", "\\2", REPLICATE))
pc_var <- pc$values %>% pull(Broken_stick) %>% (function(x) round(x*100, 1))

g_pcoa <- ggplot(pc_tab) + 
  geom_point(aes(x=PCO1, y=PCO2, colour = INDIVIDUAL), 
             size = 1.5, alpha = 0.7, shape = 16) +
  xlab(paste0("Principle co-ordinate 1 (", pc_var[1], "%)")) +
  ylab(paste0("Principle co-ordinate 2 (", pc_var[2], "%)")) +
  scale_colour_manual(values = pilot_palette, name = "Individual") + 
  theme_base

cat("done.\n")

#==============================================================================
# Make grid plot 
#==============================================================================

cat("\nAssembling output grid...")

# Extract legend
legend_out <- get_legend(g_pcoa)

# Make plot row
fig <- plot_grid(treeplot + theme(legend.position = "none"),
                 g_pcoa + theme(legend.position = "none"),
                 nrow = 1, ncol = 2, labels = "AUTO",
                 label_size = fontsize_base * fontscale_label,
                 label_fontfamily = titlefont, label_colour = "black",
                 label_fontface = "plain", vjust=1.1)

fig_out <- plot_grid(fig, legend_out, nrow = 2, rel_heights = c(1, 0.1))
                          
cat("done.\n")

#==============================================================================
# SAVE PLOT
#==============================================================================

cat("\nSaving output...")
save_fig(out_path, fig_out, out_height, out_width)
cat("...done.\n")