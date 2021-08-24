##########################################################
## AUXILIARY R FUNCTIONS FOR IGH LOCUS MAPPING WORKFLOW ##
##########################################################

#==============================================================================
# Loading packages
#==============================================================================

load_packages <- function(package_list){
    for (p in package_list){
        suppressMessages(suppressWarnings(library(p, character.only = TRUE)))
    }
}

packages_default <- c("dplyr", "stringr", "tidyr", "readr", "Biostrings",
                      "stringi", "methods", "reshape2", "entropy", "HybridMTest",
                      "data.table")
# "ggdendro", "gridExtra", , "lazyeval", "alakazam", "shazam",
# "ggplot2",
load_packages(packages_default)

#loadNamespace("alakazam")

#==============================================================================
# Auxiliary functions
#==============================================================================

q_range <- function(min_q, max_q, step_q){
  # Compute range of q-values for diversity spectra
  q <- seq(min_q, max_q, step_q)
  q0 <- (0 %in% q)
  if (!q0) {
    q <- c(0, q)
  }
  return(list(range = q, keep_0 = q0))
}

#assign_default <- function(key, default, ref = snakemake@params){
#    if (key %in% names(ref)) return(ref[[key]])
#    return(default)
#}

assign_default <- function(key, default = NULL, ref = snakemake@params){
    if (! key %in% names(ref)) return(default)
    if (ref[[key]] %in% c("None", "NULL")) return(NULL)
    return(ref[[key]])
}
