################################################################################
# Perform permutation-based Kruskal-Wallis test on bootstrapped diversity data #
################################################################################

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

# Configure output paths
output_path <- snakemake@output[[1]]

# Configure script parameters
n_permutes <- snakemake@params[["n_permutes"]]
verbose <- snakemake@params[["verbose"]]
group_within <- snakemake@params[["group_within"]]
group_between <- snakemake@params[["group_between"]]

cat("done.\n")

#===============================================================================
# Auxiliary functions
#===============================================================================

# Adapted from HybridMTest package
row.kruskal.wallis <- function (Y, grplbl) 
{
      ngenes <- dim(Y)[1]
  ugrps <- unique(grplbl)
    ngrps <- length(ugrps)
    T.mtx <- matrix(NA, ngenes, ngrps)
      n <- rep(NA, ngrps)
      tR <- apply(Y, 1, rank)
        Yrank <- t(tR)
        for (i in 1:ngrps) {
                grp.mtch <- (grplbl == ugrps[i])
            n[i] <- sum(grp.mtch)
                T.mtx[, i] <- rowMeans(Yrank[, grp.mtch])
              }
          N <- sum(n)
          k <- 12/(N * (N + 1))
            H.stat <- k * (T.mtx - (N + 1)/2)^2 %*% n
            pval <- 1 - pchisq(H.stat, ngrps - 1)
              gren.res <- grenander.ebp(unlist(pval))
              res <- cbind.data.frame(stat = H.stat, pval = pval, ebp = gren.res$ebp)
                return(res)
}

# Adapted from HybridMTest package
grenander.ebp  <- function (p) 
{
      na <- is.na(p)
  p.edf <- ecdf(p[!na])
    gren.res <- grenander(p.edf)
    if (length(gren.res$x.knots) > 1){
            gren.pdf <- approx(gren.res$x.knots, gren.res$F.knots, xout = p)$y
      } else if (length(gren.res$x.knots) == 1) {
              gren.pdf <- gren.res$F.knots[1]
        }
      gren.ebp <- min(gren.res$F.knots)/gren.pdf
      ebp.null <- pval.pdf <- rep(NA, length(p))
        ebp.null[!na] <- gren.ebp
        pval.pdf[!na] <- gren.pdf
          return(cbind.data.frame(pval.pdf = pval.pdf, ebp.null = ebp.null))
}

timetaken <- function (started.at){
    # (copied from data.table package)
    if (!inherits(started.at, "proc_time")) 
        stop("Use started.at=proc.time() not Sys.time() (POSIXt and slow)")
    format = function(secs) {
        if (secs > 60) {
            secs <- as.integer(secs)
            sprintf("%02d:%02d:%02d", secs%/%3600L, (secs%/%60L) %% 60L, secs%%60L)
        } else {
            sprintf(if (secs >= 10) "%.1fs" else "%.3fs", secs)
        }
    }
    tt = proc.time() - started.at
    paste0(format(tt[3L]), " elapsed (", format(tt[1L]), " cpu)")
}


make_permute_tabs <- function(div_tab, n_permutes, group_within,
                              group_between, verbose = FALSE){
    # Run multiple permuted KW tests and return outputs
    mclapply(1:n_permutes, function(n)
             kruskal_permute(div_tab, n, group_within, group_between, verbose),
             mc.cores = snakemake@threads) %>% bind_rows
}

kruskal_permute <- function(div_tab, permute_id, group_within,
                            group_between, verbose = FALSE){
    # Permute a diversity spectrum dataset and perform KW test
    if (verbose){
        start <- proc.time()
        msg_in <- paste0("\n\tInitialising permutation ", permute_id, ": ", date())
        cat(msg_in)
    }
    permute_tab <- permute(div_tab, permute_id, group_within, group_between)
    kw_tab <- kruskal_multi(permute_tab, group_within, group_between) %>% 
        dplyr::group_by(Q,PERMUTE) %>% summarise(P = min(P))
    if (verbose){
        msg_out <- paste0("\n\tFinished permutation ", permute_id, ": ", date(),
                          " (", timetaken(start), ")")
        cat(msg_out)
    }
    return(kw_tab)
}

permute <- function(div_tab, permute_id, group_within, group_between){
    # Permute outer-group memberships of inner groups
    group_within_old <- div_tab %>% dplyr::group_by_at(c(group_within, group_between)) %>%
        summarise %>% ungroup
    group_within_new <- group_within_old %>% 
        dplyr::mutate(GW_NEW = sample(.data[[group_within]]))
    div_new <- div_tab %>% inner_join(group_within_new, by=c(group_within, group_between)) %>%
        dplyr::select(-!!group_within) %>% dplyr::rename(!!group_within := GW_NEW) %>%
        dplyr::mutate(PERMUTE = permute_id)
    return(div_new)
}

kruskal_multi <- function(div_tab, group_within, group_between){
    # Perform parallel Kruskal-Wallis tests for each Q-value and bootstrap
    # iteration, for a given permutation of the data
    tab_spread <- div_tab %>% dplyr::arrange_at(group_between) %>%
        dplyr::select(Q, PERMUTE, ITER, !!group_between, D) %>%
        spread(!!group_between, D, sep="_") %>%
        dplyr::arrange(Q, PERMUTE, ITER)
    tab_groups <- div_tab %>% dplyr::group_by_at(c(group_between, group_within)) %>%
        summarise %>% dplyr::arrange_at(group_between) %>% pull(!!group_within)
    tab_kw <- row.kruskal.wallis(tab_spread %>% dplyr::select(-Q, -ITER, -PERMUTE),
                                tab_groups) %>% as_tibble %>%
        dplyr::mutate(Q = tab_spread$Q, PERMUTE=tab_spread$PERMUTE,
               ITER = tab_spread$ITER)
    #cat(colnames(tab_kw))
    tab_kw <- dplyr::rename(tab_kw, P = "pval", H = "stat", EBP = "ebp")
    return(tab_kw)
}

#===============================================================================
# Import data
#===============================================================================

cat("\nImporting bootstrap table...")
col <- cols(ITER = "i", N_GROUP = "i", Q = "d", D = "d",
            .default = "c")
diversities_bootstrapped <- suppressMessages(read_tsv(input_path, col_types = col))
cat("done.\n")
cat("Input dimensions:", dim(diversities_bootstrapped), "\n")

#===============================================================================
# Perform analysis
#===============================================================================

# Perform kruskal-wallis tests for real data
cat("\nComputing KW p-values from true data...")
tab_pkw_true  <- diversities_bootstrapped %>% mutate(PERMUTE = 0) %>%
    kruskal_multi(group_within, group_between)
cat("done.\n")

# Perform kruskal-wallis tests for real data
cat("\nPreparing comparison traces...")
cat("\n\tMinimum p-value...")
tab_pkw_min   <- tab_pkw_true %>% dplyr::group_by(Q) %>% summarise(P = min(P)) %>%
    dplyr::mutate(METHOD = "KW (minimum P-value)")
cat("\n\tAverage diversity...")
tab_pkw_avg_d <- diversities_bootstrapped %>% mutate(PERMUTE = 0) %>%
    dplyr::group_by_at(c("Q", "PERMUTE", group_within, group_between)) %>%
    summarise(D = mean(D)) %>% dplyr::mutate(ITER = 0) %>% ungroup %>%
    kruskal_multi(group_within, group_between) %>% dplyr::mutate(METHOD = "KW (mean D-value)")
cat("\n\tAggregating comparison traces...")
tab_pkw_true_all <- bind_rows(tab_pkw_min, tab_pkw_avg_d) %>%
    dplyr::select(-H, -EBP, -PERMUTE, -ITER)
cat("\n...done.\n")

# Run permutations and get KW p-values for each
cat("\nNumber of permutations to run:", n_permutes, "\n")
cat("\nExecuting permutation tests...")
tab_permutes <- make_permute_tabs(diversities_bootstrapped, n_permutes,
                                  group_within, group_between, verbose)
cat("...permutation tests complete.\n")

# Compute permutation p-values for real data
cat("\nComputing real p-values from permutation data...")
tab_pkw_permute_raw <- tab_pkw_min %>% dplyr::rename(P_TRUE = P) %>%
    dplyr::select(-METHOD) %>%
    inner_join(tab_permutes, by="Q") %>% dplyr::mutate(P_LESS = P <= P_TRUE)
tab_pkw_permute <- tab_pkw_permute_raw %>%
    dplyr::group_by(Q) %>% summarise(P = mean(P_LESS)) %>% 
    dplyr::mutate(METHOD = "KW permutation test")
cat("\ndone.\n")

# Concatenate different methods
tab_pkw_out <- bind_rows(tab_pkw_true_all, tab_pkw_permute)

#===============================================================================
# Write output
#===============================================================================

cat("Writing p-values to file...")
write_tsv(tab_pkw_out, output_path)
cat("done.\n")
