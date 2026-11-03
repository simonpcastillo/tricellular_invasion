# =============================================================================
# PVN vs GLOBAL (WHOLE-TISSUE) ENRICHMENT
# =============================================================================
# "Global" here = ALL cells in the tissue, NOT just cells inside EC
# neighbourhoods. This is a different reference from "vs other PVNs", and needs
# a different null:
#
#   vs other PVNs  -> shuffle PVN labels among ECs        (pvn_enrichment_tests.R)
#   vs global      -> shuffle CELL labels across tissue   (this file)
#
# The vs-global null preserves the whole-tissue composition and each EC's
# neighbour count/positions, but randomises WHICH cells carry which label. It
# asks: does the perivascular space concentrate this cell type / tumour state
# beyond what a random spatial arrangement of the same tissue would give?
#
# Marginal (strata = Image) vs region-adjusted (strata = c("Image","region")):
# the marginal null breaks ALL spatial structure, so it flags anything that is
# spatially clustered at all (incl. purely regional gradients). The
# region-adjusted null permutes WITHIN region, so it isolates enrichment in the
# immediate perivascular space beyond the regional background -- that is the
# interpretable "perivascular vs tissue" test. Report both.
#
# Source AFTER section 6 of 01_pvn_discovery.R (needs df0 with class/Classification
# /cx/cy/gid/Image/region, and cnt with Image/e/PVN).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(spatstat.geom)
})

# ---- EC-anchored non-EC neighbour pairs, restricted to the retained EC set --
build_pairs <- function(df0, cnt, radius) {
  retained <- unique(cnt[, .(Image, e)])
  pl <- rbindlist(lapply(unique(df0$Image), function(img) {
    s  <- df0[Image == img]
    pp <- ppp(s$cx, s$cy, window = owin(range(s$cx), range(s$cy)), marks = s$class)
    cp <- closepairs(pp, r = radius, twice = FALSE)
    if (!length(cp$i)) return(NULL)
    I <- s$gid[cp$i]; J <- s$gid[cp$j]; ec <- s[class == "Endothelial", gid]
    keep <- I %in% ec | J %in% ec; I <- I[keep]; J <- J[keep]
    data.table(Image = img,
               e  = fifelse(I %in% ec, I, J),
               nb = fifelse(I %in% ec, J, I))
  }))
  pl <- pl[df0$class[nb] != "Endothelial"]     # non-EC neighbours only
  pl <- pl[retained, on = .(Image, e), nomatch = 0]
  pl[]
}

.strat_vec <- function(df0, strata) {
  if (length(strata) == 1) as.character(df0[[strata]])
  else do.call(paste, c(lapply(strata, function(s) as.character(df0[[s]])), sep = "||"))
}

# permute a per-cell label vector within stratum (only over the given cell set)
.permute_within <- function(lab, strat, cell_idx) {
  out <- lab
  s <- split(cell_idx, strat[cell_idx])
  for (ix in s) out[ix] <- lab[sample(ix)]
  out
}

# =============================================================================
# (ii-global) MICROGLIA: perivascular fraction vs whole-tissue fraction
# =============================================================================
# Per EC: microglia neighbours / all non-EC neighbours (same denominator as
# prop[,MG]). Niche statistic = mean of that over the niche's ECs. Null shuffles
# the microglia label among ALL non-EC cells within stratum.
micro_vs_global <- function(df0, cnt, pl, mg_class,
                            strata = c("Image"), n_perm = 2000, seed = 1) {
  ekey  <- paste(pl$Image, pl$e)
  ef    <- factor(ekey); lv <- levels(ef)
  tot_e <- as.numeric(table(ef))                       # neighbours per EC
  pvn_e <- cnt[match(lv, paste(Image, e)), PVN]
  is_mg <- df0$class == mg_class
  nb_pos <- pl$nb                                      # gid of neighbour cell
  strat  <- .strat_vec(df0, strata)
  nonec  <- which(df0$class != "Endothelial")

  niche_means <- function(mgvec) {
    mg_e <- rowsum(as.numeric(mgvec[nb_pos]), ef)[, 1]
    frac <- mg_e / tot_e
    tapply(frac, pvn_e, mean)
  }
  obs <- niche_means(is_mg)
  set.seed(seed)
  null <- replicate(n_perm, niche_means(.permute_within(is_mg, strat, nonec)))
  .summarise(obs, null, n_perm, id = "PVN",
             extra = data.table(mg_frac_obs = as.numeric(obs),
                                global_frac = mean(is_mg[nonec])))
}

# =============================================================================
# (iii-global) NEFTEL: perivascular within-tumour state fraction vs whole-tumour
# =============================================================================
# Per EC: state neighbours / all tumour neighbours. Null shuffles the tumour
# STATE label among ALL tumour cells within stratum (tumour identity fixed, so
# the reference is the whole-tumour state composition).
neftel_vs_global <- function(df0, cnt, pl, states,
                             strata = c("Image"), n_perm = 2000, seed = 1) {
  ncls  <- df0$Classification[pl$nb]
  plt   <- pl[ncls %in% states]                        # tumour neighbours only
  ncls  <- ncls[ncls %in% states]
  ekey  <- paste(plt$Image, plt$e)
  ef    <- factor(ekey); lv <- levels(ef)
  tot_e <- as.numeric(table(ef))                       # tumour neighbours per EC
  pvn_e <- cnt[match(lv, paste(Image, e)), PVN]
  nb_pos <- plt$nb
  strat  <- .strat_vec(df0, strata)
  tum    <- which(df0$Classification %in% states)
  state_lab <- df0$Classification

  niche_means <- function(lab) {
    labv <- lab[nb_pos]
    setNames(lapply(states, function(st) {
      c_e  <- rowsum(as.numeric(labv == st), ef)[, 1]
      tapply(c_e / tot_e, pvn_e, mean)
    }), states)
  }
  obs <- niche_means(state_lab)
  set.seed(seed)
  null_raw <- replicate(n_perm,
                        niche_means(.permute_within(state_lab, strat, tum)),
                        simplify = FALSE)
  gf <- sapply(states, function(st) mean(df0$Classification[tum] == st))
  rbindlist(lapply(states, function(st) {
    nn <- sapply(null_raw, function(z) z[[st]])          # niches x n_perm
    if (is.null(dim(nn))) nn <- matrix(nn, nrow = length(obs[[st]]))
    r <- .summarise(obs[[st]], nn, n_perm, id = "PVN",
                    extra = data.table(state = st,
                                       state_frac_obs = as.numeric(obs[[st]]),
                                       global_frac = gf[[st]]))
    r
  }))[order(state, -effect)]
}

# ---- shared summariser: two-sided p vs the permutation expectation ---------
.summarise <- function(obs, null, n_perm, id = "PVN", extra = NULL) {
  nms <- names(obs); obs <- as.numeric(obs)
  if (is.null(dim(null))) null <- matrix(null, nrow = length(obs))
  nm  <- as.numeric(rowMeans(null, na.rm = TRUE))      # null expectation per niche
  eff <- obs - nm
  p   <- vapply(seq_along(obs), function(k)
    (1 + sum(abs(null[k, ] - nm[k]) >= abs(obs[k] - nm[k]) - 1e-12, na.rm = TRUE)) /
      (n_perm + 1), numeric(1))
  out <- data.table(id = nms, n_null_mean = nm,
                    effect = eff, p_perm = p)
  setnames(out, "id", id)
  if (!is.null(extra)) out <- cbind(out, extra)
  out[, p_BH := p.adjust(p_perm, "BH")]
  out[]
}
