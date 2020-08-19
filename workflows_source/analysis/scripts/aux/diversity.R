# Custom diversity functions based on those from Alakazam

loadNamespace("alakazam")

#------------------------------------------------------------------------------
# Private/auxiliary functions copied from Alakazam
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

q_range <- function(min_q, max_q, step_q){
  # Compute range of q-values for diversity spectra
  q <- seq(min_q, max_q, step_q)
  q0 <- (0 %in% q)
  if (!q0) {
    q <- c(0, q)
  }
  return(list(range = q, keep_0 = q0))
}

ci_scale <- function(ci){
  # Compute z-score and SD scale from confidence interval, assuming gaussian
  ci_z <- ci + (1 - ci)/2
  return(qnorm(ci_z))
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
                                    clone_field){
  # For each inner group in an outer group, compute bootstrap counts
  # and return a concatenated counts table
  bootstraps <- tibble()
  gtab <- group_tab[group_tab[[group_within_name]] == group_within_value,]
  for (g in gtab %>% pull(group_between_name) %>% unique){
    sample_tab <- bootstrap_abundance(nboot, clone_tab, group_tab, group_within_name, 
                                      group_within_value, group_between_name,
                                      g, clone_field)
    bootstraps <- bind_rows(bootstraps, sample_tab)
  }
  return(bootstraps)
}

compute_outer_bootstrap <- function(nboot, clone_tab, group_tab, group_within_name,
                                    group_between_name, clone_field){
  # Compute inner bootstrap values for each outer group and collate values
  # in a single table
  bootstraps <- tibble()
  for (G in group_tab %>% pull(group_within_name) %>% unique){
    sample_tab <- compute_inner_bootstrap(nboot, clone_tab, group_tab, group_within_name,
                                          G, group_between_name, clone_field)
    bootstraps <- bind_rows(bootstraps, sample_tab)
  }
  return(bootstraps)
}

#------------------------------------------------------------------------------
# New diversity functions
#------------------------------------------------------------------------------

compute_solo_diversity_bootstrap <- function(bootstraps, group_within,
                                             group_between, clone_field, Q){
    # Compute solo diversity spectrum within each inner group and bootstrap rep
    bs_tab <- bootstraps %>%
        dplyr::group_by(!!sym(group_within), !!sym(group_between), ITER)
    div_tab_split <- tibble()
    for (q in Q){
        dt <- dplyr::summarise(bs_tab, Q = q, D = calcDiversity(N, q = q),
                        N_GROUP = 1)
        div_tab_split <- bind_rows(div_tab_split, dt)
    }
    return(div_tab_split)
}

compute_gamma_diversity_bootstrap <- function(bootstraps, group_within,
                                              clone_field, Q){
    # Compute gamma diversity within each outer group and bootstrap rep
    bs_tab <- bootstraps %>% 
      dplyr::group_by(!!sym(group_within), !!sym(group_between), ITER) %>%
      dplyr::summarise(N = sum(N)) %>%
      dplyr::group_by(!!sym(group_within),  ITER)
    div_tab <- tibble()
    for (q in Q){
        dt <- dplyr::summarise(bs_tab, Q = q, D = calcDiversity(N, q = q))
        div_tab <- bind_rows(div_tab, dt)
    }
    return(div_tab)
}

compute_alpha_diversity_bootstrap <- function(bootstraps, group_within,
                                              group_between, clone_field, Q){
    # Compute alpha diversity within each outer group and bootstrap rep
    bs_tab <- bootstraps %>%
      dplyr::group_by(!!sym(group_within), !!sym(group_between), ITER)
    div_tab_split <- tibble()
    for (q in Q){
        dt <- dplyr::summarise(bs_tab, Q = q, D = calcDiversity(N, q = q))
        div_tab_split <- bind_rows(div_tab_split, dt)
    }
    div_tab <- div_tab_split %>% 
        dplyr::group_by(!!sym(group_within), ITER, Q) %>% 
        dplyr::summarise(D = ifelse(dplyr::first(Q) != 1,
                             mean(D^(1-dplyr::first(Q)))^(1/(1-dplyr::first(Q))),
                             exp(mean(log(D)))),
                  N_GROUP = n())
    return(div_tab)
}

compute_beta_diversity_bootstrap <- function(gamma_div_tab, alpha_div_tab,
                                             group_within){
    # Compute beta diversity spectrum within each outer group and bootstrap rep,
    # from pre-computed gamma and alpha spectra
    full_join(gamma_div_tab, alpha_div_tab, by = c(group_within, "ITER", "Q"),
              suffix = c("_GAMMA", "_ALPHA")) %>%
        dplyr::mutate(D= D_GAMMA/D_ALPHA) %>% dplyr::select(-D_GAMMA, -D_ALPHA)
}

summarise_bootstrap_diversity <- function(div_tab, group_within, ci, keep_0,
                                          group_between = NA,
                                          empirical_bounds = TRUE){
  #' dplyr::summarise diversity table across bootstrap repeats
  # Calculate empirical quantiles
  quantile_upper <- 1-(1-ci)/2
  quantile_lower <- (1-ci)/2
  # Compute summary table
  summ_tab <- div_tab %>% dplyr::group_by(!!sym(group_within), Q) %>%
    {if ("N_GROUP" %in% colnames(.)) dplyr::group_by(., N_GROUP,
                                              add = TRUE) else . } %>%
    {if (!is.na(group_between)) dplyr::group_by(., !!sym(group_between),
                                         add = TRUE) else . } %>%
    dplyr::summarise(D_SD = sd(D),
              D_LOWER_EMPIRICAL = quantile(D, quantile_lower),
              D_UPPER_EMPIRICAL = quantile(D, quantile_upper),
              D = mean(D)
              ) %>%
    dplyr::mutate(D_LOWER_NORMAL = D - ci_scale(ci) * D_SD,
           D_UPPER_NORMAL = D + ci_scale(ci) * D_SD)
    # Add evenness information
    n = ifelse(!is.na(group_between), 2, 1)
    summ_tab_evenness <- summ_tab %>%
        dplyr::group_by(!!sym(group_within)) %>%
        {if (!is.na(group_between)) dplyr::group_by(., !!sym(group_between),
                                              add = TRUE) else . } %>%
        filter(Q==0) %>% transmute(D_0 = D) %>%
        full_join(summ_tab,
                  by = c(group_within, group_between)[1:n]) %>%
        dplyr::mutate(E = D/D_0, 
               E_LOWER_EMPIRICAL = D_LOWER_EMPIRICAL/D_0, 
               E_UPPER_EMPIRICAL = D_UPPER_EMPIRICAL/D_0,
               E_LOWER_NORMAL = D_LOWER_NORMAL/D_0, 
               E_UPPER_NORMAL = D_UPPER_NORMAL/D_0
               ) %>%
        dplyr::select(-D_0)
    # Drop 0-order row if not in desired Q-range
    if (keep_0) return(summ_tab_evenness)
    return(sum_tab_evenness %>% filter(Q != 0))
}

#------------------------------------------------------------------------------
# Wrapper functions for computing all three grouped diversities
#------------------------------------------------------------------------------

compute_diversity_bootstraps <- function(tab, group_within, group_between,
                                         clone_field = "CLONE", nboot = 2000,
                                         copy_field = NULL, min_n = 30,
                                         max_n = NULL, uniform = TRUE){
    # Compute bootstraps for diversity-curve estimation
    # Configure test variables
    data = check_data(tab, group_within, group_between,
                      clone_field, copy_field)
    # Compute clone copy numbers for each group based on row numbers 
    # (if copy column is NULL) or copy column values
    clone_tab <- make_clone_tab(data, group_within, group_between,
                                clone_field, copy_field)
    # Count sequences in each group
    group_tab <- make_group_tab(clone_tab, group_within, group_between, min_n,
                                max_n, uniform)
    # Generate table of bootstrap counts over all groupings
    bootstraps <- compute_outer_bootstrap(nboot, clone_tab, group_tab,
                                          group_within, group_between,
                                          clone_field)
    return(bootstraps)
}

compute_diversity_grouped <- function(bootstraps, group_within, group_between,
                                      clone_field = "CLONE", min_q = 0, 
                                      max_q = 4, step_q = 0.05, ci = 0.95){
    # Compute grouped alpha/beta/gamma diversity spectra from bootstraps
    # Compute Q value range
    q <- q_range(min_q, max_q, step_q)
    # Compute gamma and alpha diversity for each bootstrap replicate
    gamma_div_bs <- compute_gamma_diversity_bootstrap(bootstraps, group_within,
                                                      clone_field, q$range)
    alpha_div_bs <- compute_alpha_diversity_bootstrap(bootstraps, group_within,
                                                      group_between, 
                                                      clone_field, q$range)
    # Compute beta diversity from gammma and alpha
    beta_div_bs <- compute_beta_diversity_bootstrap(gamma_div_bs,
                                                    alpha_div_bs, group_within)
    # Compute averages and distributions over all bootstrap replicates
    alpha_div <- summarise_bootstrap_diversity(alpha_div_bs, group_within, ci,
                                               q$keep_0, NA)
    beta_div <- summarise_bootstrap_diversity(beta_div_bs, group_within, ci,
                                              q$keep_0, NA)
    gamma_div <- summarise_bootstrap_diversity(gamma_div_bs, group_within, ci,
                                               q$keep_0, NA)
    # Compute pairwise, Q-wise p-values for differences in diversity
    alpha_test <- compute_diversity_pvalues(alpha_div_bs, group_within)
    gamma_test <- compute_diversity_pvalues(gamma_div_bs, group_within)
    beta_test  <- compute_diversity_pvalues(beta_div_bs, group_within)
    # Collapse all tables together (one dplyr::summarised, one for p-values)
    all_div <- bind_rows(alpha_div %>% dplyr::mutate(DIVTYPE = "alpha"),
                         gamma_div %>% dplyr::mutate(DIVTYPE = "gamma"),
                         beta_div %>% dplyr::mutate(DIVTYPE = "beta"))
    all_test <- bind_rows(alpha_test %>% dplyr::mutate(DIVTYPE = "alpha"),
                          gamma_test %>% dplyr::mutate(DIVTYPE = "gamma"),
                          beta_test %>% dplyr::mutate(DIVTYPE = "beta"))
    return(list(div=all_div, test=all_test))
}


compute_diversity_solo <- function(bootstraps, group_within, group_between,
                                   clone_field = "CLONE", min_q = 0, 
                                   max_q = 4, step_q = 0.05, ci = 0.95){
    # Compute grouped alpha/beta/gamma diversity spectra from bootstraps
    # Compute Q value range
    q <- q_range(min_q, max_q, step_q)
    # Compute solo diversity for each bootstrap replicate
    solo_div_bs <- compute_solo_diversity_bootstrap(bootstraps, group_within,
                                                    group_between, 
                                                    clone_field, q$range)
    # Compute averages and distributions over all bootstrap replicates
    solo_div <- summarise_bootstrap_diversity(solo_div_bs, group_within, ci,
                                              q$keep_0, group_between) %>%
        dplyr::mutate(DIVTYPE = "solo")
    solo_test <- compute_diversity_pvalues(solo_div_bs, group_within) %>%
      dplyr::mutate(DIVTYPE = "solo")
    return(list(div=solo_div, test=solo_test))
}

compute_diversity_pvalues <- function(bootstrap_diversities, group_within){
  #' Compute the pairwise statistical significance of differences
  #' in bootstrapped diversity values, based on ECDF of
  #' bootstrap diversities (derived from Alakazam helperTest function).
  group_pairs <- combn(unique(bootstrap_diversities[[group_within]]), 2)
  qvals <- sort(unique(bootstrap_diversities[["Q"]]))
  pvalue_tab <- tibble()
  # Iterate over groups and q-values
  for (n in 1:ncol(group_pairs)){
    group_pair <- group_pairs[,n]
    pair_list <- list()
    for (q in qvals){
      mat1 <- bootstrap_diversities %>%
        filter(!!sym(group_within) == group_pair[1], Q == q) %>%
        pull(D)
      mat2 <- bootstrap_diversities %>%
        dplyr::filter(!!sym(group_within) == group_pair[2], Q == q) %>%
        pull(D)
      if (mean(mat1) >= mean(mat2)) { 
        g_delta <- mat1 - mat2 
      } else { 
        g_delta <- mat2 - mat1 
      }
      p <- ecdf(g_delta)(0)
      p <- ifelse(p <= 0.5, p * 2, (1 - p) * 2)
      pair_list[[as.character(q)]] <- list(DELTA_MEAN=mean(g_delta), 
                                           DELTA_SD=sd(g_delta), 
                                           PVALUE=p)
    }
    pvalue_tab <- bind_rows(pvalue_tab,
                            bind_rows(pair_list, .id = "Q") %>% 
                              dplyr::mutate(GROUP1 = group_pair[1], GROUP2 = group_pair[2])) %>%
      dplyr::select(GROUP1, GROUP2, everything())
  }
  return(pvalue_tab)
}
