#######################################################################
# Collapse clones into consensus sequence for SHM modelling          ##
#######################################################################

#==============================================================================
# Auxiliary functions
#==============================================================================

load_packages <- function(package_list){
  for (p in package_list){
    suppressMessages(suppressWarnings(library(p, character.only = TRUE)))
  }
}

assign_default <- function(key, default = NULL, ref = snakemake@params){
  if (! key %in% names(ref)) return(default)
  if (ref[[key]] %in% c("None", "NULL")) return(NULL)
  return(ref[[key]])
}

#==============================================================================
# Loading packages
#==============================================================================

packages_default <- c("dplyr", "readr", "shazam")
# "alakazam", "stringr", "tidyr", "reshape2", "entropy", "stringi", "methods", "Biostrings",
load_packages(packages_default)

#==============================================================================
# Preamble
#==============================================================================

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Configure input paths
input_path <- snakemake@input[[1]]

# Configure output paths
output_path <- snakemake@output[[1]]

# Configure script parameters
clone_field <- assign_default("clone_field", "CLONE")
sequence_field <- assign_default("sequence_field", "SEQUENCE_IMGT")
germline_field <- assign_default("germline_field", "GERMLINE_IMGT_D_MASK")
method <- assign_default("method", "mostCommon")
include_ambiguous <- assign_default("include_ambiguous", FALSE)
break_ties_stochastic <- assign_default("break_ties_stochastic", FALSE)
break_ties_columns <- assign_default("break_ties_stochastic", NULL)
expanded_db <- assign_default("expanded_db", FALSE)
exclude_na <- assign_default("exclude_na", FALSE)
mu_freq_field <- assign_default("mu_freq_field", NULL)
min_freq <- assign_default("min_freq", NULL)
region_def_id <- assign_default("region_def", NULL)
if (is.null(region_def_id)){
    region_def <- NULL
} else {
    region_def <- get(region_def_id)
}

cat("done.\n")

#==============================================================================
# Import data
#==============================================================================

cat("\nImporting Change-O table...")
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(input_path, col_types = col))
cat("done.\n")
cat("Input dimensions:", dim(tab), "\n")

#==============================================================================
# Handle NA sequences
#==============================================================================

nna <- sum(is.na(tab[[clone_field]]))
cat(ifelse(exclude_na, "\nRemoving", "\nRenaming"), nna,
    "sequences without clonal IDs...")
if (exclude_na){
    tab <- tab[!is.na(tab[[clone_field]]),]
} else {
    tab[is.na(tab[[clone_field]]),][[clone_field]] <- paste(NA, seq(nna), sep = "_")
}
cat("done.\n")
exp_out <- tab %>% group_by_(.dots = clone_field) %>% summarise %>% nrow
cat("Expected output dimensions:", exp_out, ncol(tab) + 3, "\n")

#==============================================================================
# Compute clone counts and compute consensus sequences
#==============================================================================

# Compute clone counts
cat("\nComputing clonal counts...")
tab_cc <- tab %>% group_by_(.dots = clone_field) %>%
  mutate(CLNCOUNT = n()) %>% ungroup()
cat("done.\n")

# Compute consensus sequences across clones
cat("\nCollapsing sequence clones...")
clncounts <- tab_cc %>% pull("CLNCOUNT") %>% unique %>% sort
tabs_col <- lapply(clncounts, function(c)
                   collapseClones(tab_cc %>% filter(CLNCOUNT == c),
                                  cloneColumn = clone_field,
                                  sequenceColumn = sequence_field,
                                  germlineColumn = germline_field,
                                  muFreqColumn = mu_freq_field,
                                  regionDefinition = region_def,
                                  method = method,
                                  minimumFrequency = min_freq,
                                  includeAmbiguous = include_ambiguous,
                                  breakTiesStochastic = break_ties_stochastic,
                                  breakTiesByColumns = break_ties_columns,
                                  expandedDb = expanded_db,
                                  nproc = snakemake@threads))
tab_collapsed <- bind_rows(tabs_col) %>%
    rename(CLONAL_SEQUENCE = clonal_sequence,
           CLONAL_GERMLINE = clonal_germline)
cat("done.\n")

# Remove residual ambiguous characters
if (!include_ambiguous){
    cat("\nRemoving residual ambiguous characters...")
    chars_allowed <- c("A","T","G","C","N",".","-")
    regex_disallowed <- paste0("[^", paste(chars_allowed, collapse = ""), "]")
    tab_collapsed <- tab_collapsed %>%
        mutate(CLONAL_SEQUENCE = gsub(regex_disallowed, "N", CLONAL_SEQUENCE),
               CLONAL_GERMLINE = gsub(regex_disallowed, "N", CLONAL_GERMLINE))
    cat("done.\n")
}

#==============================================================================
# Save output
#==============================================================================

cat("\nOutput dimensions:", dim(tab_collapsed))
cat("\nWriting output table to file...")
write_tsv(tab_collapsed, output_path)
cat("done.\n")
