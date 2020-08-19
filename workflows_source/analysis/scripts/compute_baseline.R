######################################################################
# Count selection pressure from SHM frequencies with BASELINe        #
######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Specify defaults
getnull <- function(varname){
    !exists(varname) || get(varname) %in% c("None", "NULL")
    }
if (!exists("sequence_field")) sequence_field <- "CLONAL_SEQUENCE"
if (!exists("germline_field")) germline_field <- "CLONAL_GERMLINE"
if (!exists("test_statistic")) test_statistic <- "focused"
if (getnull("mutation_def")) mutation_def <- NULL
if (!exists("calc_stats")) calc_stats <- FALSE
if (getnull("region_def")) {
    region_def <- IMGT_V
} else {
    region_def <- get(region_def)
}

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(MU_COUNT_CDR_R = "d", MU_COUNT_CDR_S = "d",
            MU_COUNT_FWR_R = "d", MU_COUNT_FWR_S = "d",
            MU_EXPECTED_CDR_R = "d", MU_EXPECTED_CDR_S = "d",
            MU_EXPECTED_FWR_R = "d", MU_EXPECTED_FWR_S = "d",
            .default = col_character())
tab <- suppressMessages(read_tsv(inpath_tab, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Read SHM model
write_log("\nImporting mutation model...", newline = FALSE)
mod <- readRDS(inpath_mod)
log_done()

# Compute selection
write_log("\nComputing selection pressures...", newline = FALSE)
baseline <- calcBaseline(tab,
                         sequenceColumn = sequence_field,
                         germlineColumn = germline_field,
                         testStatistic  = test_statistic,
                         regionDefinition = region_def,
                         targetingModel = mod, 
                         mutationDefinition = mutation_def,
                         calcStats = calc_stats,
                         nproc = snakemake@threads)

log_done()

# Write output
write_log("\nSaving to output...", newline = FALSE)
saveRDS(baseline, outpath)
log_done()
