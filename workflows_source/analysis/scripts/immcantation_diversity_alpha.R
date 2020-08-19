############################################################
# Compute grouped alpha diversities from solo diversities ##
############################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")
if (!exists("ci")) ci <- 0.95 # Use standard default confidence intervals
if (!exists("nboot")) nboot <- 2000 # Use standard bootstrap # repeats
if (!exists("uniform")) uniform <- TRUE # Downsample counts by default
if (!exists("max_n")) max_n <- NULL # Downsample counts by default
if (copy == "NULL") copy <- NULL
if (max_n %in% c("None", "NULL")) max_n <- NULL

# Read Change-O table
write_log("\nImporting Change-O table...", newline = FALSE)
tab <- suppressMessages(readChangeoDb(inpath))
log_done()
write_log("Input dimensions:", dim(tab))

# Generate diversity S4 object
tab <- tab[!is.na(tab[[group_within]]),]
tab <- tab[!is.na(tab[[group_between]]),]
tab <- tab[!is.na(tab[[clone_field]]),]
write_log("Computing alpha-spectra over",
          paste0(group_within, "/", group_between, "..."), newline = FALSE)
div <- rarefyDiversityAlpha2(tab = tab, group_within = group_within,
                            group_between = group_between, clone_field = clone_field,
                            copy_field = copy, uniform = uniform,
                            min_q = min_q, max_q = max_q, step_q = step_q,
                            min_n = min_n, max_n = max_n, ci = ci,
                            nboot = nboot)
log_done()

# Write output
write_log("Writing abundance object to file...")
dput(div, outpath)
log_done()

# Wildcards: group_within, group_between, set
# Params: min_q, max_q, step_q, min_n, max_n, gamma, clone_field, aux, uniform,
#         ci, nboot, ...
# Write output
