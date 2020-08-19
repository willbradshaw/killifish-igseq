######################################################################
# Generate an IGOR J-anchor file from an IgBLAST auxiliary file      #
######################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)
write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Import sequences
write_log("\nReading IMGT-gapped V-sequences...", newline = FALSE)
v_imgt <- readDNAStringSet(inpath)
log_done()
write_log(length(v_imgt), "sequences read.")

# Infer anchors
write_log("\nInferring anchor positions of ungapped sequences...",
          newline = FALSE)
anchors <- v_imgt %>% subseq(1,209) %>% as.character %>% 
  gsub("\\.", "", .) %>% str_count
log_done()

# Create anchor table
write_log("\nCreating anchor table...", newline = FALSE)
anchor_tab <- tibble(NAME = names(v_imgt), ANCHOR = anchors)
log_done()

# Write output
write_log("\nOutput dimensions:", dim(anchor_tab))
write_log("\nWriting output table to file...", newline = FALSE)
write_delim(anchor_tab, outpath, delim=";", col_names = FALSE)
log_done()
