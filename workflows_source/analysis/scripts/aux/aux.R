##########################################################
## AUXILIARY R FUNCTIONS FOR IGH LOCUS MAPPING WORKFLOW ##
##########################################################

#-----------------------------------------------------------------------------
# Loading packages (also for these functions)
#-----------------------------------------------------------------------------

load_packages <- function(package_list){
    for (p in package_list){
        suppressMessages(suppressWarnings(library(p, character.only = TRUE)))
    }
}

packages_default <- c(#"genbankr", "BSgenome", "rentrez", "GenomicRanges", 
                      "dplyr", "stringr", "tidyr", "readr", "Biostrings",
                      "alakazam", "shazam", "stringi", "methods", "ggplot2",
                      "gridExtra", "ggdendro", "reshape2", "entropy",
                      "lazyeval", "HybridMTest")
load_packages(packages_default)

#-----------------------------------------------------------------------------
# Install and load missing packages from miniCRAN
#-----------------------------------------------------------------------------

if (! is.null(snakemake@params[["minicran_path"]])){
    packages_install <- snakemake@params[["minicran_install"]]
    repo <- snakemake@params[["minicran_path"]] %>% normalizePath() %>%
        paste0("file://", .)
    if (! is.null(snakemake@params[["minicran_save"]])){
        if (!dir.exists(snakemake@params[["minicran_save"]])){
            dir.create(snakemake@params[["minicran_save"]],
                       recursive = TRUE)
        }
        lib = snakemake@params[["minicran_save"]]
        suppressMessages(suppressWarnings(install.packages(packages_install,
                                                           repos = repo,
                                                           lib = lib)))
    } else {
        suppressMessages(suppressWarnings(install.packages(packages_install,
                                                           repos=repo)))
    }
    load_packages(packages_install)
}

#-----------------------------------------------------------------------------
# Source other auxiliary packages
#-----------------------------------------------------------------------------

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "dplyr.R"))
source(file.path(aux_dir, "diversity.R"))
# TODO: Source other aux packages for igseq if needed

#-----------------------------------------------------------------------------
# Parsing Snakemake object
#-----------------------------------------------------------------------------

abspath <- function(path, dir=getwd()){
    # Add a directory to a relative file path to make it absolute
    return(file.path(dir, path))
}

assign_global <- function(name, value, abs=FALSE){
    # Assign a value to a name in the global environment
    # HANDLE WITH CARE! Not for normal functions - just snakemake assignment
    if (abs) value <- abspath(value)
    assign(name, value, envir = .GlobalEnv)
}

parse_named <- function(object, prefix, abs=FALSE){
    # Parse elements of an object with named elements into separate objects
    for (n in names(object)){
        name <- ifelse(nchar(prefix) > 0, paste0(prefix, "_", n), n)
        assign_global(name, object[[n]], abs)
    }
}

parse_unnamed <- function(object, prefix, abs=FALSE){
    # Parse elements of an object without named elements into separate objects
    if (length(object) == 1){
        assign_global(prefix, object[[1]], abs)
    } else {
        assign_global(paste0(prefix, "s"), object, abs)
    }
}

parse_snakemake_slot <- function(snakemake, slot_name, prefix, abs=FALSE){
    # Parse elements of a snakemake object slot into separate objects
    o <- slot(snakemake, slot_name)
    if (length(o) == 0) return()
    if (is.null(names(o))) return(parse_unnamed(o, prefix, abs))
    o <- o[str_count(names(o)) > 0]
    return(parse_named(o, prefix, abs))
}

parse_snakemake <- function(snakemake){
    # Parse slots from global snakemake object into separate objects
    parse_snakemake_slot(snakemake, "wildcards", "", FALSE)
    parse_snakemake_slot(snakemake, "params", "", FALSE)
    parse_snakemake_slot(snakemake, "input", "inpath", FALSE)
    parse_snakemake_slot(snakemake, "output", "outpath", TRUE)
    parse_snakemake_slot(snakemake, "log", "logpath", FALSE)
}

# TODO: Add parse reporting (e.g. log("Input path (GenBank file):", inpath_gb))

#-----------------------------------------------------------------------------
# Log script progress
#-----------------------------------------------------------------------------

clear_log <- function(path=get("logpath")){
    # Delete a pre-existing log file
    if (file.exists(path)) file.remove(path)
}

write_log <- function(..., newline=TRUE, path=get("logpath")){
    # Append messages to log file
    cat(c(...), ifelse(newline, "\n", ""), file=path, append=TRUE)
}

log_done <- function(newline=FALSE, path=get("logpath")) {
    write_log("done.", ifelse(newline, "\n", ""), path=path)
}

#-----------------------------------------------------------------------------
# Fixed-pattern grep
#-----------------------------------------------------------------------------

fgrep <- function(...) grep(..., fixed=TRUE)
fgrepl <- function(...) grepl(..., fixed=TRUE)
