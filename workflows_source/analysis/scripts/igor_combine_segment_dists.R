######################################################################
# Melt and combine segment distributions into a single table        ##
######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Read input tables
write_log("\nReading segment distribution tables...", newline = FALSE)
col <- cols(.default = col_character())
tabs <- lapply(inpaths, function(p)
               suppressMessages(read_tsv(p, col_types = col)))
log_done()

# Melt and combine tables
write_log("\nMelting and concatenating tables...", newline = FALSE)
tabs_melted <- lapply(tabs, function(t)
    t %>% melt(id.vars = c("id", "p"), variable.name = "segment_types",
               value.name = "segments") %>%
        mutate(segment_types = toupper(segment_types))
        )
tab_out <- bind_rows(tabs_melted)
log_done()

# Write output
write_log("\nOutput dimensions:", dim(tab_out))
write_log("\nWriting output table to file...", newline = FALSE)
write_tsv(tab_out, outpath)
log_done()
