# =============================================================================
# PVN ENRICHMENT STATISTICS  (add-on to 01_pvn_discovery.R)
# =============================================================================
# Answers three questions with EC-level permutation tests (the EC neighbourhood
# is the unit; region is available as a stratum). Source AFTER section 9 of the
# main script so cnt / prop / nef / df0 exist.
#
#   (ii)  microglia_niche_enrichment()  -- does a PVN have more microglia than
#                                          the other PVNs / than the global mean
#   (iii) neftel_niche_enrichment()     -- each state in each PVN vs other PVNs
#                                          (== vs global; see note)
#         neftel_group_contrast()       -- microglia-enriched PVNs vs the rest
#
# All tests are TWO-SIDED and permute the grouping label WITHIN a stratum
# (default Image; pass strata=c("Image","region") for the region-adjusted
# version). BH-correct across the family you care about.
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

# ---- assemble one per-EC table ---------------------------------------------
# Returns: Image, e, PVN, region, mg_frac, neftel_<state>, tum_tot
build_ec_table <- function(cnt, prop, nef, df0, states = NEFTEL, mg_class = NULL) {
  ec <- data.table(Image = cnt$Image, e = cnt$e, PVN = cnt$PVN,
                   region = df0$region[cnt$e])
  if (is.null(mg_class)) {
    cand <- grep("microglia", colnames(prop), ignore.case = TRUE, value = TRUE)
    mg_class <- if (length(cand)) cand[1] else NA_character_
  }
  ec[, mg_frac := if (!is.na(mg_class)) prop[, mg_class] else NA_real_]
  attr(ec, "mg_class") <- mg_class
  nf <- nef[, c("Image", "e", paste0("neftel_", states), "tum_tot"), with = FALSE]
  ec <- merge(ec, nf, by = c("Image", "e"), all.x = TRUE)
  ec[]
}

.stratum <- function(ec, strata) {
  if (length(strata) == 1) as.character(ec[[strata]])
  else do.call(paste, c(lapply(strata, function(s) as.character(ec[[s]])), sep = "||"))
}

# ---- engine: each level of `label` vs the rest, within-stratum permutation --
# nl (per-level count) is preserved by within-stratum shuffling, so the
# vs-global effect = mean_level - grand_mean is a fixed monotone rescale of the
# vs-rest effect: the permutation p-value is IDENTICAL for the two references.
# We report both effect sizes and the single shared p.
.perm_vs_rest <- function(value, label, stratum, n_perm = 2000, seed = 1) {
  ok <- is.finite(value) & !is.na(label) & !is.na(stratum)
  value <- value[ok]; label <- as.character(label)[ok]; stratum <- stratum[ok]
  n <- length(value); grand <- mean(value); levs <- sort(unique(label))
  ml <- tapply(value, label, mean)[levs]
  nl <- as.numeric(table(factor(label, levs)))
  eff_rest <- as.numeric(ml - (grand * n - ml * nl) / (n - nl))
  eff_glob <- as.numeric(ml - grand)
  sidx <- split(seq_len(n), stratum)
  ge <- numeric(length(levs)); set.seed(seed)
  for (b in seq_len(n_perm)) {
    perm <- label
    for (ix in sidx) perm[ix] <- sample(label[ix])
    mlp  <- tapply(value, perm, mean)[levs]
    effp <- as.numeric(mlp - (grand * n - mlp * nl) / (n - nl))
    ge   <- ge + (abs(effp) >= abs(eff_rest) - 1e-12)
  }
  data.table(level = levs, n_ec = nl, mean = as.numeric(ml), grand = grand,
             effect_vs_rest = eff_rest, effect_vs_global = eff_glob,
             p_perm = (1 + ge) / (n_perm + 1))
}

# ---- engine: two-group contrast (enriched vs not), within-stratum permutation
.perm_group <- function(value, grp, stratum, n_perm = 2000, seed = 1) {
  ok <- is.finite(value) & !is.na(grp) & !is.na(stratum)
  value <- value[ok]; grp <- as.logical(grp)[ok]; stratum <- stratum[ok]
  obs <- mean(value[grp]) - mean(value[!grp])
  sidx <- split(seq_along(value), stratum); set.seed(seed); ge <- 0
  for (b in seq_len(n_perm)) {
    g <- grp; for (ix in sidx) g[ix] <- sample(grp[ix])
    ge <- ge + (abs(mean(value[g]) - mean(value[!g])) >= abs(obs) - 1e-12)
  }
  list(effect = obs, p_perm = (1 + ge) / (n_perm + 1),
       mean_grp = mean(value[grp]), mean_rest = mean(value[!grp]),
       n_grp = sum(grp), n_rest = sum(!grp))
}

# =============================================================================
# (ii) MICROGLIA: is a PVN enriched vs the other PVNs / vs the global mean?
# =============================================================================
microglia_niche_enrichment <- function(ec, strata = c("Image"),
                                       n_perm = 2000, seed = 1) {
  stopifnot("mg_frac" %in% names(ec))
  d <- ec[!is.na(PVN) & is.finite(mg_frac)]
  r <- .perm_vs_rest(d$mg_frac, d$PVN, .stratum(d, strata), n_perm, seed)
  setnames(r, "level", "PVN")
  r[, p_BH := p.adjust(p_perm, "BH")]
  setnames(r, c("mean", "grand"), c("mg_frac_mean", "mg_frac_global"))
  r[order(-effect_vs_rest)]
}

# parametric cross-check: binomial GLM on microglia COUNTS (respects the
# differing neighbour totals per EC). One-vs-rest Wald test per PVN.
microglia_niche_glm <- function(cnt, mg_class, covars = c("Image")) {
  d <- copy(cnt)[!is.na(PVN)]
  stopifnot(mg_class %in% names(d))
  d[, mg := get(mg_class)]
  levs <- sort(unique(d$PVN))
  out <- rbindlist(lapply(levs, function(k) {
    d[, inPVN := factor(PVN == k, c(FALSE, TRUE))]
    f  <- reformulate(c("inPVN", covars), response = "cbind(mg, total - mg)")
    fit <- tryCatch(glm(f, data = d, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) return(data.table(PVN = k, log_OR = NA_real_, p = NA_real_))
    co <- summary(fit)$coefficients; row <- grep("inPVNTRUE", rownames(co))
    data.table(PVN = k, log_OR = co[row, 1], p = co[row, 4])
  }))
  out[, p_BH := p.adjust(p, "BH")][order(-log_OR)]
}

# =============================================================================
# (iii-a) NEFTEL state enrichment: each state in each PVN vs other PVNs
#         (== vs global; shared p — see .perm_vs_rest note)
# =============================================================================
neftel_niche_enrichment <- function(ec, states = NEFTEL, strata = c("Image"),
                                    n_perm = 2000, seed = 1) {
  res <- rbindlist(lapply(states, function(st) {
    v <- ec[[paste0("neftel_", st)]]
    d <- data.table(v = v, PVN = ec$PVN, strat = .stratum(ec, strata))
    d <- d[!is.na(PVN) & is.finite(v)]
    r <- .perm_vs_rest(d$v, d$PVN, d$strat, n_perm, seed)
    setnames(r, "level", "PVN"); r[, state := st][]
  }))
  res[, p_BH := p.adjust(p_perm, "BH")]          # BH across all PVN x state
  setcolorder(res, c("state", "PVN", "n_ec", "mean", "grand",
                     "effect_vs_rest", "effect_vs_global", "p_perm", "p_BH"))
  res[order(state, -effect_vs_rest)]
}

# =============================================================================
# (iii-b) NEFTEL: microglia-ENRICHED PVNs vs the rest, per state
# =============================================================================
neftel_group_contrast <- function(ec, enriched_pvns, states = NEFTEL,
                                  strata = c("Image"), n_perm = 2000, seed = 1) {
  stopifnot(length(enriched_pvns) >= 1)
  res <- rbindlist(lapply(states, function(st) {
    v   <- ec[[paste0("neftel_", st)]]
    grp <- ec$PVN %in% enriched_pvns
    g   <- .perm_group(v, grp, .stratum(ec, strata), n_perm, seed)
    data.table(state = st, mean_enriched = g$mean_grp, mean_other = g$mean_rest,
               effect = g$effect, n_enriched = g$n_grp, n_other = g$n_rest,
               p_perm = g$p_perm)
  }))
  res[, p_BH := p.adjust(p_perm, "BH")][order(-abs(effect))]
}

# ---- convenience: run a test marginal (Image) AND region-adjusted, merged ----
run_dual <- function(fun, ..., key) {
  m <- fun(..., strata = c("Image"))
  a <- fun(..., strata = c("Image", "region"))
  num <- setdiff(intersect(names(m), names(a)), key)
  merge(m[, c(key, num), with = FALSE], a[, c(key, num), with = FALSE],
        by = key, suffixes = c("_marg", "_adj"))
}
