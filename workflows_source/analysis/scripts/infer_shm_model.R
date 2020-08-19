#######################################################################
# Infer an SHM targeting model from a Change-O sequence DB            #
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
if (!exists("model")) model <- "S"
if (!exists("sequence_field")) sequence_field <- "CLONAL_SEQUENCE"
if (!exists("germline_field")) germline_field <- "CLONAL_GERMLINE"
if (!exists("v_call_field")) v_call_field <- "BEST_V_CALL"
if (!exists("mult_mut")) mult_mut <- "independent"
if (!exists("min_mut")) min_mut <- 50
if (!exists("min_seq_mut")) min_seq_mut <- 500
if (!exists("model_name")) model_name <- ""
if (!exists("model_desc")) model_desc <- ""
if (!exists("model_species")) model_species <- ""
if (!exists("model_cite")) model_cite <- ""
if (getnull("model_date")) model_date <- NULL

write_log("Model type:", model)

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
col <- cols(.default = col_character())
tab <- suppressMessages(read_tsv(inpath, col_types = col))
log_done()
write_log("Input dimensions:", dim(tab))

# Compute targeting model
write_log("\nComputing targeting model...", newline = FALSE)
tarmod <- createTargetingModel(tab, 
                               model = model,
                               sequenceColumn = sequence_field,
                               germlineColumn = germline_field,
                               vCallColumn = v_call_field,
                               multipleMutation = mult_mut,
                               minNumMutations = min_mut,
                               minNumSeqMutations = min_seq_mut,
                               modelName = model_name,
                               modelDescription = model_desc,
                               modelSpecies = model_species,
                               modelCitation = model_cite,
                               modelDate = model_date)
log_done()

# Combine into output table and write to output
write_log("Writing output to file...")
saveRDS(tarmod, outpath)
log_done()
