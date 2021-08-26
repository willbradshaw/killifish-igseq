#######################################################################
# Compute bootstrap replicates from a Change-O table, saving each     #
# replicate separately.                                               #
#######################################################################

#==============================================================================
# Preamble
#==============================================================================

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Source packages and auxiliary functions
source("scripts/aux.R")

# Configure input paths
input_path <- snakemake@input[[1]]

# Configure output directory (from checkpoint rule)
output_dir <- snakemake@output[[1]]

# Configure script parameters
min_n <- snakemake@params[["min_n"]]
max_n <- snakemake@params[["max_n"]]
nboot <- snakemake@params[["nboot"]]
clone_field <- snakemake@params[["clone_field"]]
group_within <- snakemake@params[["group_within"]]
group_between <- snakemake@params[["group_between"]]
out_pattern <- snakemake@params[["out_pattern"]]
replicate_field <- snakemake@params[["replicate_field"]]
verbose <- TRUE # TODO: Add to params
if (max_n %in% c("None", "NULL")) max_n <- NULL

# Configure output file paths
output_path_base <- file.path(output_dir, out_pattern)
process_path <- function(n){
    pattern <- paste0("\\{", replicate_field, "\\}")
    return(gsub(pattern, n, output_path_base))
}

cat("done.\n")

#==============================================================================
# Auxiliary functions
#==============================================================================

#------------------------------------------------------------------------------
# Private/auxiliary functions copied from Alakazam (v. 0.2.10)
#------------------------------------------------------------------------------

checkColumns <- function (data, columns, logic = c("all", "any")){
    logic <- match.arg(logic)
    data_names <- names(data)
    if (logic == "all") {
        for (f in columns) {
            if (!(f %in% data_names)) { 
                msg <- paste("The column", f, "was not found")
                return(msg)
            }
        }
        for (f in columns) {
            if (all(is.na(data[[f]]))) {
                msg <- paste("The column", f, "contains no data")
                return(msg)
            }
        }
    } else if (logic == "any") {
        if (!any(columns %in% data_names)) {
            msg <- paste("Input must contain at least one of the columns:", 
                         paste(columns, collapse = ", "))
            return(msg)
        }
        columns_found <- columns[columns %in% data_names]
        invalid <- sapply(columns_found, function(f) all(is.na(data[[f]])))
        if (all(invalid)) {
            msg <- paste("None of the columns", paste(columns_found, 
                                                      collapse = ", "), "contain data")
            return(msg)
        }
    }
    return(TRUE)
}

inferUnseenCount <- function(x) {
    x <- x[x >= 1]
    n <- sum(x)
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    if (f2 > 0) {
        f0 <- ceiling(((n - 1) * f1^2) / (n * 2 * f2))
    } else {
        f0 <- ceiling(((n - 1) * f1 * (f1 - 1)) / (n * 2))
    }
    return(f0)
}

inferUnseenAbundance <- function(x) {
    x <- x[x >= 1]
    # Coverage
    rC1 <- calcCoverage(x, r=1)
    # Unseen count
    f0 <- inferUnseenCount(x)
    # Assign unseen relative abundance
    p <- rep((1 - rC1) / f0, f0)
    return(p)
}

adjustObservedAbundance <- function(x) {
    x <- x[x >= 1]
    n <- sum(x)
    # Coverage
    rC1 <- calcCoverage(x, r=1)
    # Calculate tuning parameter
    lambda <- (1 - rC1) / sum(x/n * exp(-x))
    # Define adjusted relative abundance
    p <- x/n * (1 -  lambda * exp(-x))
    return(p)
}

calcCoverage <- function (x, r = 1){
  if (r == 1) {
    return(calcChao1Coverage(x))
  }
  x <- x[x >= 1]
  n <- sum(x)
  fr <- sum(x == r)
  fs <- sum(x == r + 1)
  if (fr == 0) {
    stop("Cannot calculate coverage of order ", r, ". No abundance data with count=", 
         r, ".")
  }
  if (fs == 0) {
    stop("Cannot calculate coverage of order ", r, ". No abundance data with count=", 
         r + 1, ".")
  }
  a <- factorial(r) * fr/sum(x[x >= r]^r)
  b <- ((n - r) * fr/((n - r) * fr + (r + 1) * fs))^r
  rC <- 1 - a * b
  return(rC)
}

calcChao1Coverage <- function(x) {
  x <- x[x >= 1]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  
  if (f2 > 0) {
    rC1 <- 1 - (f1 / n) * (((n - 1) * f1) / ((n - 1) * f1 + 2 * f2))
  } else {
    rC1 <- 1 - (f1 / n) * (((n - 1) * (f1 - 1)) / ((n - 1) * (f1 - 1) + 2))
  }
  
  return(rC1)
}

#------------------------------------------------------------------------------
# New auxiliary functions
#------------------------------------------------------------------------------

check_data <- function(data, group_within, group_between, clone_field,
                       copy_field){
  # Filter and validate input data for diversity computation
  # Confirm input data is a data frame
  if (!is.data.frame(data)) {
    stop("Input data is not a data.frame")
  }
  # Check that group, clone and copy columns in data are valid
  check <- checkColumns(data, c(group_within, group_between, clone_field, 
                                copy_field))
  if (check != TRUE) { stop(check)  }
  data <- data[!is.na(data[[clone_field]]),]
  data <- data[!is.na(data[[group_within]]),]
  data <- data[!is.na(data[[group_between]]),]
  # Reformat group and clone columns as strings
  data[[group_within]] <- as.character(data[[group_within]])
  data[[group_between]] <- as.character(data[[group_between]])
  data[[clone_field]] <- as.character(data[[clone_field]])
  return(data)
}

make_clone_tab <- function(data, group_within, group_between, clone_field, copy_field){
  # Compute copy numbers for each clone group, annotated by sample group and 
  # normalised by copy field if specified
  if (is.null(copy_field)) {
    clone_tab <- dplyr::group_by(data, !!sym(group_within), !!sym(group_between),
                          !!sym(clone_field)) %>%
      dplyr::summarise(COUNT = n())
  } else {
    clone_tab <- dplyr::group_by(data, !!sym(group_within), !!sym(group_between),
                          !!sym(clone_field)) %>%
      dplyr::summarise(COUNT = sum(!!sym(copy_field)))
  }
  return(clone_tab)
}

make_group_tab <- function(clone_tab, group_within, group_between, min_n,
                           max_n, uniform){
    # dplyr::summarise and process clone tab by group membership to get group counts
    group_tab <- clone_tab %>% dplyr::group_by(!!sym(group_within), !!sym(group_between)) %>%
      dplyr::summarize(SEQUENCES = sum(COUNT, na.rm = TRUE)) %>%
      dplyr::mutate(KEEP = SEQUENCES >= min_n)
    # Warn if some groups excluded due to sequence number
    if (any(!group_tab$KEEP)){
        group_tab_out <- group_tab %>% filter(!KEEP)
        group_str_out <- paste(group_tab_out[[group_within]], 
                               group_tab_out[[group_between]], sep = "//")
        warning("Not all groups passed threshold min_n=", min_n, ". ",
                "Excluded: ",  paste(group_str_out, collapse = ", "), ".")
    }
    # Exclude too-small groups and compute sample sizes
    group_tab <- group_tab %>% filter(KEEP) %>% dplyr::select(-KEEP) %>% ungroup()
    if (uniform) {
        group_tab <- dplyr::mutate(group_tab, NSAM = min(SEQUENCES, max_n))
    } else if (is.null(max_n)){
        group_tab <- dplyr::mutate(group_tab, NSAM = SEQUENCES)
    } else {
        group_tab <- dplyr::mutate(group_tab, NSAM = pmin(SEQUENCES, max_n))
    }
    return(group_tab)
}

bootstrap_abundance <- function(nboot, clone_tab, group_tab, group_within_name, 
                                group_within_value, group_between_name,
                                group_between_value, clone_field){
  # Generate independent bootstrap samples for a given group combination
  # in a clone table
  gtab <- group_tab[group_tab[[group_within_name]] == group_within_value &
                      group_tab[[group_between_name]] == group_between_value,]
  ctab <- clone_tab[clone_tab[[group_between_name]] == group_between_value &
                      clone_tab[[group_within_name]] == group_within_value,]
  # Get sequence number for resampling
  n <- gtab$NSAM
  # Get abundances of each clone in group
  abund_obs <- ctab$COUNT
  # Infer abundances of observed and unseen clones (Chao method)
  p1 <- adjustObservedAbundance(abund_obs)
  p2 <- inferUnseenAbundance(abund_obs)
  abund_inf <- c(p1, p2)
  # Specify clone names
  n1 <- ctab[[clone_field]]
  n2 <- paste("UNKNOWN", group_between_value, seq_along(p2), sep = "_")
  names_inf <- c(n1, n2)
  # Use inferred clone distribution and resampling size to 
  # generate independent bootstrap samples
  sample_mat <- rmultinom(nboot, n, abund_inf)
  sample_tab <- melt(sample_mat, varnames = c(clone_field, "ITER"),
                       value.name = "N")
  sample_tab[[clone_field]] <- names_inf[sample_tab[[clone_field]]]
  sample_tab[[group_within_name]] <- group_within_value
  sample_tab[[group_between_name]] <- group_between_value
  return(sample_tab)
}

compute_inner_bootstrap <- function(nboot, clone_tab, group_tab, group_within_name,
                                    group_within_value, group_between_name,
                                    clone_field, verbose){
  # For each inner group in an outer group, compute bootstrap counts
  # and return a concatenated counts table
  bootstraps <- tibble()
  gtab <- group_tab[group_tab[[group_within_name]] == group_within_value,]
  gvec <- gtab %>% pull(group_between_name) %>% unique %>% sample
  for (g in gvec){
    if (verbose){
        start <- proc.time()
        msg_in <- paste0("\n\tComputing bootstraps for ", group_within_value,
                         "//", g, "...")
        cat(msg_in)
    }
    sample_tab <- bootstrap_abundance(nboot, clone_tab, group_tab, group_within_name,
                                      group_within_value, group_between_name,
                                      g, clone_field)
    bootstraps <- bind_rows(bootstraps, sample_tab)
    if (verbose){
        cat("done. (", timetaken(start), ")")
    }
  }
  gc(verbose=FALSE)
  return(bootstraps)
}

compute_outer_bootstrap <- function(nboot, clone_tab, group_tab, group_within_name,
                                    group_between_name, clone_field, verbose){
  # Compute inner bootstrap values for each outer group and collate values
  # in a single table
  bootstraps <- tibble()
  for (G in group_tab %>% pull(group_within_name) %>% unique){
    sample_tab <- compute_inner_bootstrap(nboot, clone_tab, group_tab, group_within_name,
                                          G, group_between_name, clone_field, verbose)
    bootstraps <- bind_rows(bootstraps, sample_tab)
  }
  gc(verbose=FALSE)
  return(bootstraps)
}

compute_diversity_bootstraps <- function(tab, group_within, group_between,
                                         clone_field = "CLONE", nboot = 2000,
                                         copy_field = NULL, min_n = 30,
                                         max_n = NULL, uniform = TRUE,
                                         verbose = TRUE){
    # Compute bootstraps for diversity-curve estimation
    # Configure test variables
    if (verbose) cat("\n\tChecking data...")
    data = check_data(tab, group_within, group_between,
                      clone_field, copy_field)
    if (verbose) cat("done.")
    # Compute clone copy numbers for each group based on row numbers 
    # (if copy column is NULL) or copy column values
    if (verbose) cat("\n\tGenerating clone table...")
    clone_tab <- make_clone_tab(data, group_within, group_between,
                                clone_field, copy_field)
    if (verbose) cat("done.")
    # Count sequences in each group
    if (verbose) cat("\n\tGenerating group table...")
    group_tab <- make_group_tab(clone_tab, group_within, group_between, min_n,
                                max_n, uniform)
    if (verbose) cat("done.")
    # Generate table of bootstrap counts over all groupings
    bootstraps <- compute_outer_bootstrap(nboot, clone_tab, group_tab,
                                          group_within, group_between,
                                          clone_field, verbose)
    gc(verbose=FALSE)
    return(bootstraps)
}

#==============================================================================
# Import data
#==============================================================================

cat("\nImporting Change-O table...")
col <- cols(DUPCOUNT = "n", CONSCOUNT = "n",
            .default = col_character())
tab <- suppressMessages(read_tsv(input_path, col_types = col))
cat("done.\n")
cat("Input dimensions:", dim(tab), "\n")

#==============================================================================
# Generate and save bootstrap replicates
#==============================================================================

# Make output directory
if (!dir.exists(output_dir)){
    cat("\nCreating output directory...")
    dir.create(output_dir, recursive=TRUE)
    cat("done.\n")
}

# Run bootstraps
cat("\nComputing bootstrap table...")
for (n in 1:nboot){
    # Compute replicate
    cat("\nComputing bootstrap replicate", n, "...")
    bootstraps <- compute_diversity_bootstraps(tab, group_within, group_between,
                                               clone_field, 1, NULL,
                                               min_n, max_n, TRUE, verbose)
    bootstraps_out <- mutate(bootstraps, ITER = n)
    cat("\n...done.")
    cat("\nOutput dimensions:", dim(bootstraps))
    # Save replicate
    out_path_rep <- process_path(n)
    cat("\nOutput path:", out_path_rep)
    cat("\nWriting bootstrap replicate", n, "to file...")
    write_tsv(bootstraps_out, out_path_rep)
    cat("\n...done.\n")
}
