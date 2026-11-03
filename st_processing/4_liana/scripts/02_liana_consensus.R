################################################################################
# 02_liana_consensus.R
# Project: Perivascular invasion paper
# Author: Simon Castillo
# Last update: 19/09/2026
################################################################################
##
##

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(purrr)
})

# Resolve paths independently of the current working directory.
.args_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(.args_file)) stop("Run this script with Rscript")
.script_path <- normalizePath(sub("^--file=", "", .args_file[[1]]), mustWork = TRUE)
PROJECT_ROOT <- dirname(dirname(.script_path))
.config_file <- Sys.getenv("LIANA_CONFIG", file.path(PROJECT_ROOT, "config", "config.R"))
if (!file.exists(.config_file)) stop("Configuration not found: ", .config_file)
source(.config_file, local = FALSE)

ROOT <- LIANA_RESULTS
REPS <- unname(SAMPLES)
modes <- MODES
n_rep <- length(REPS)
if (n_rep < 2L) stop("Consensus requires at least two replicate samples")
KEY <- c("source", "target", "ligand.complex", "receptor.complex")

for(mode in modes){ 
  OUT_DIR <- file.path(ROOT, "consensus", mode)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # ---- 1. read all replicates ----------------------------------------------
  long <- bind_rows(lapply(REPS, function(s)
    read_csv(file.path(ROOT, s, mode ,"liana_aggregate.csv"), 
             show_col_types = FALSE) |>
      mutate(replicate = s)))
  
  if (SHARED_ONLY) {
    shared <- Reduce(intersect, lapply(split(long, long$replicate),
                                       function(d) union(d$source, d$target)))
    dropped <- setdiff(union(long$source, long$target), shared)
    message("Groups in all replicates: ", paste(sort(shared), collapse = ", "))
    if (length(dropped))
      message("Dropped (absent from >=1 replicate): ", 
              paste(sort(dropped), collapse = ", "))
    long <- filter(long, source %in% shared, target %in% shared)
  }
  
  if (anyDuplicated(long[, c(KEY, "replicate")]))
    stop("Duplicate LR pair × replicate rows in mode ", mode)
  long <- long |>
    group_by(replicate) |>
    mutate(rank_pct = percent_rank(aggregate_rank),
           p = pmax(aggregate_rank, P_FLOOR)) |>
    ungroup()
  
  # ---- 2. per-interaction summary across replicates ------------------------
  # Fisher combination below reproduces the historical calculation.
  # aggregate_rank is a rank-aggregation score, NOT necessarily a calibrated
  # p-value. Treat p_fisher/q_fisher as heuristic until independently validated.
  wide <- long |>
    group_by(across(all_of(KEY))) |>
    summarise(n_tested      = n(),
              p_max         = max(p),
              fisher_chi    = -2 * sum(log(p)),
              mean_rank_pct = mean(rank_pct),
              mean_LRscore  = if ("sca.LRscore" %in% names(long))
                mean(sca.LRscore, na.rm = TRUE) else NA_real_,
              mean_spec     = if ("natmi.edge_specificity" %in% names(long))
                mean(natmi.edge_specificity, na.rm = TRUE) else NA_real_,
              .groups = "drop") |>
          mutate(
            p_fisher = pchisq(
              fisher_chi,
              df = 2 * n_tested,
              lower.tail = FALSE
            ),
            q_fisher = p.adjust(p_fisher, method = "BH")
          )
  
  # per-replicate columns alongside, for inspection
  wide <- long |>
    select(all_of(KEY), 
           replicate, 
           aggregate_rank, 
           rank_pct,
           any_of(c("sca.LRscore", "cellphonedb.pvalue"))) |>
    pivot_wider(id_cols = all_of(KEY), 
                names_from = replicate,
                values_from = -all_of(c(KEY, "replicate")), 
                names_sep = "__") |>
    right_join(wide, by = KEY)
  
  # The original script calculated Fisher's p/q values but then overwrote
  # that filter with hard-coded CellPhoneDB p-values for 7000 and 7003.
  # Make the historical choice explicit and support a q-based alternative.
  cp_cols <- grep("^cellphonedb\\.pvalue__", names(wide), value = TRUE)
  has_cpdb <- if (length(cp_cols)) {
    cp <- as.matrix(wide[, cp_cols, drop = FALSE])
    rowSums(!is.na(cp) & cp < CPDB_P) > 0L
  } else rep(FALSE, nrow(wide))
  if (CONSENSUS_FILTER == "legacy_cpdb_any") {
    if (!length(cp_cols)) stop("No CellPhoneDB p-value columns for historical filter")
    consensus <- wide[has_cpdb, , drop = FALSE]
  } else if (CONSENSUS_FILTER == "fisher_q") {
    consensus <- wide[is.finite(wide$q_fisher) & wide$q_fisher <= FDR, , drop = FALSE]
  } else if (CONSENSUS_FILTER == "fisher_q_and_cpdb_any") {
    consensus <- wide[is.finite(wide$q_fisher) & wide$q_fisher <= FDR & has_cpdb, , drop = FALSE]
  } else stop("Unknown CONSENSUS_FILTER: ", CONSENSUS_FILTER)
  if (!is.na(CONSENSUS_MIN_LRSCORE)) {
    consensus <- consensus[is.finite(consensus$mean_LRscore) &
                           consensus$mean_LRscore >= CONSENSUS_MIN_LRSCORE, , drop = FALSE]
  }
  consensus <- arrange(consensus, q_fisher, desc(mean_spec))
  message("Consensus mode=", CONSENSUS_FILTER, ", retained ", nrow(consensus), " pairs")

  wide <- wide |> as.data.frame()
  consensus <- consensus |> as.data.frame()
  
  write_csv(wide ,      file.path(OUT_DIR, "liana_replicate_wide.csv"))
  write_csv(consensus, file.path(OUT_DIR, "liana_consensus.csv"))
  
  cat(sprintf("testable in all %d replicates: %d of %d interactions\n",
              n_rep, sum(wide$n_tested == n_rep), nrow(wide)))
  cat(sprintf("consensus set: %d interactions across %d source->target pairs\n",
              nrow(consensus), n_distinct(paste(consensus$source, 
                                                consensus$target))))
  print(consensus |> 
          count(source, target, name = "n") |> 
          arrange(desc(n)) |> 
          head(10))
  
  # ---- 3. replicate concordance QC -----------------------------------------
  pairs_qc <- combn(REPS, 2, simplify = FALSE) |>
    map_dfr(function(pr) {
      a <- filter(long, replicate == pr[1]) 
      b <- filter(long, replicate == pr[2])
      m <- inner_join(a, b, by = KEY, suffix = c(".a", ".b"))
      siga <- with(m, aggregate_rank.a <= 0.01)
      sigb <- with(m, aggregate_rank.b <= 0.01)
      tibble(rep_a = pr[1], rep_b = pr[2],
             n_common = nrow(m),
             spearman_rho = cor(m$rank_pct.a, m$rank_pct.b, method = "spearman"),
             n_sig_a = sum(siga), n_sig_b = sum(sigb),
             jaccard_sig = sum(siga & sigb) / max(sum(siga | sigb), 1))
    })
  
  pairs_qc <- pairs_qc |> as.data.frame()
  write_csv(pairs_qc, file.path(OUT_DIR, "liana_qc_concordance.csv"))
  print(pairs_qc)
  
  # ---- 4. figures -----------------------------------------------------------
  pr <- REPS[1:2]
  cons_key <- consensus |> 
    select(all_of(KEY)) |> 
    mutate(consensus = TRUE)
  
  sc <- long |>
    filter(replicate %in% pr) |>
    select(all_of(KEY), replicate, aggregate_rank) |>
    pivot_wider(id_cols = all_of(KEY), names_from = replicate,
                values_from = aggregate_rank) |>
    filter(!is.na(.data[[pr[1]]]), !is.na(.data[[pr[2]]])) |>
    left_join(cons_key, by = KEY) |>
    mutate(consensus = ifelse(is.na(consensus), "not in consensus", "consensus")) |> 
    as.data.frame()
  
  p_sc <- ggplot(sc, aes(-log10(pmax(.data[[pr[1]]], P_FLOOR)),
                         -log10(pmax(.data[[pr[2]]], P_FLOOR)),
                         colour = consensus)) +
    geom_abline(slope = 1, 
                intercept = 0, 
                linetype = 2, 
                colour = "grey60") +
    geom_point(alpha = 0.55, 
               size = 0.9) +
    scale_colour_manual(values = c(consensus = "#B2182B",
                                   `not in consensus` = "grey70")) +
    labs(x = paste0("-log10 aggregate_rank  (", pr[1], ")"),
         y = paste0("-log10 aggregate_rank  (", pr[2], ")"), 
         colour = NULL,
         subtitle = sprintf("Spearman rho = %.3f on %d interactions tested in both",
                            pairs_qc$spearman_rho[1], pairs_qc$n_common[1])) +
    theme_bw(base_size = 11)
  
  ggsave(file.path(OUT_DIR, "liana_concordance.pdf"), 
         plot = p_sc,
         width = 7.5, height = 5.5)
  
  top <- consensus |> slice_head(n = NTOP_DOT) |>
    mutate(interaction = paste(ligand.complex, receptor.complex, sep = " -> "),
           interaction = factor(interaction, levels = rev(unique(interaction)))) |> 
    as.data.frame()
  
  if (nrow(top) == 0L) {
    message("No consensus pairs in mode ", mode, "; dotplot skipped")
    next
  }
  p_dot <- ggplot(top, aes(target, interaction,
                           size = -log10(q_fisher), 
                           colour = mean_LRscore)) +
    geom_point() +
    facet_grid(~ source, 
               scales = "free_x", 
               space = "free_x") +
    scale_colour_viridis_c(option = "C") +
    scale_size_continuous(range = c(1.5, 6)) +
    labs(x = "target (receptor)", y = NULL,
         size = expression(-log[10]~q), 
         colour = "mean\nLRscore",
         title = sprintf("Top %d reproducible interactions (n = %d replicates)",
                         nrow(top), n_rep)) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
  ggsave(file.path(OUT_DIR, "liana_consensus_dotplot.pdf"), plot = p_dot,
         width = 5 + 1.1 * n_distinct(top$source), height = 8, limitsize = FALSE)
  
  message("Done. Outputs in ", normalizePath(OUT_DIR))
}
print('End of 02_liana_consensus.R -------------------------------------------')


