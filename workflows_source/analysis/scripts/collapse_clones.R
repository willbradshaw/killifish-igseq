#######################################################################
# Collapse clones into consensus sequence for SHM modelling          ##
#######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Specify defaults
getnull <- function(varname){
    !exists(varname) || get(varname) %in% c("None", "NULL")
    }
if (!exists("clone_field")) clone_field <- "CLONE"
if (!exists("sequence_field")) sequence_field <- "SEQUENCE_IMGT"
if (!exists("germline_field")) germline_field <- "GERMLINE_IMGT_D_MASK"
if (getnull("mu_freq_field")) mu_freq_field <- NULL
if (!exists("method")) method <- "mostCommon"
if (getnull("min_freq")) min_freq <- NULL
if (!exists("include_ambiguous")) include_ambiguous <- FALSE
if (!exists("break_ties_stochastic")) break_ties_stochastic <- FALSE
if (getnull("break_ties_columns")) break_ties_columns <- NULL
if (!exists("expanded_db")) expanded_db <- FALSE
if (!exists("exclude_na")) exclude_na <- FALSE # Retain sequneces without clonal ID?
if (getnull("region_def")) {
    region_def <- NULL
} else {
    region_def <- get(region_def)
}

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Handle NA sequences
nna <- sum(is.na(tab[[clone_field]]))
if (exclude_na){
    write_log("\nRemoving", nna,
              "sequences without clonal IDs...", newline = FALSE)
    tab <- tab[!is.na(tab[[clone_field]]),]
    log_done()
} else {
    write_log("\nRenaming", nna,
              "sequences without clonal IDs...", newline = FALSE)
    tab[is.na(tab[[clone_field]]),][[clone_field]] <- paste(NA, seq(nna), sep = "_")
    log_done()
}
exp_out <- tab %>% group_by_(.dots = clone_field) %>% summarise %>% nrow
write_log("Expected output dimensions:", exp_out, ncol(tab) + 3)

# Compute clone counts
write_log("\nComputing clonal counts...", newline = FALSE)
tab_cc <- tab %>% group_by_(.dots = clone_field) %>%
  mutate(CLNCOUNT = n()) %>% ungroup()
log_done()

# Compute consensus sequences across clones
write_log("\nCollapsing sequence clones...", newline = FALSE)
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
tab_collapsed <- bind_rows(tabs_col)
log_done()

# Remove residual ambiguous characters
if (!include_ambiguous){
    chars_allowed <- c("A","T","G","C","N",".","-")
    regex_disallowed <- paste0("[^", paste(chars_allowed, collapse = ""), "]")
    tab_collapsed <- tab_collapsed %>%
        mutate(CLONAL_SEQUENCE = gsub(regex_disallowed, "N", CLONAL_SEQUENCE),
               CLONAL_GERMLINE = gsub(regex_disallowed, "N", CLONAL_GERMLINE))
}

# Combine into output table and write to output
#tab_out <- bind_rows(tab_solo, tab_multi)
write_log("\nOutput dimensions:", dim(tab_collapsed))
write_log("Writing output to file...")
write_tsv(tab_collapsed, outpath)
log_done()
