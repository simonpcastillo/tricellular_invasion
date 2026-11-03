# =============================================================================
# Standalone bulk RNA-seq module
# Configuration: set BULKRNASEQ_CONFIG to a config.R path (see run.sh).
# Unspecified settings fall back to the defaults defined in this script.
# =============================================================================
CONFIG_PATH <- Sys.getenv("BULKRNASEQ_CONFIG", unset = "")
.config_env <- new.env(parent = globalenv())
if (nzchar(CONFIG_PATH)) {
  if (!file.exists(CONFIG_PATH)) stop("Configuration file missing: ", CONFIG_PATH)
  sys.source(CONFIG_PATH, envir = .config_env)
}
cfg <- function(name, default) {
  if (exists(name, envir = .config_env, inherits = FALSE))
    get(name, envir = .config_env, inherits = FALSE)
  else default
}

################################################################################
# Project: Perivascular invasion paper
# Author: Simon Castillo
# Last update: 17/09/2026
################################################################################
##
## REGION SIGNATURES + FINAL FIGURE  
## 0. Preamble ----

suppressPackageStartupMessages({
  library(limma); library(edgeR)
  library(ggplot2); library(tidyr)      
  library(ggnewscale); library(patchwork)
  library(dplyr)                        
})

.masked <- function() {
  v <- intersect(c("desc", "slice", "rename", "first", "count", "select",
                   "union", "intersect", "setdiff", "collapse"),
                 getNamespaceExports("dplyr"))
  v[vapply(v, function(f) { w <- find(f)
  !length(w) || w[1] != "package:dplyr" }, logical(1))]
}
if (length(.masked())) {
  detach("package:dplyr", unload = FALSE, character.only = TRUE)
  suppressPackageStartupMessages(library(dplyr))
  if (length(.bad <- .masked()))
    stop("dplyr verbs still masked after re-attaching: ",
         paste(.bad, collapse = ", "), "\n  restart R and source this file first.")
  cat("NOTE: dplyr was masked by an S4 package and has been re-attached on top\n")
}
ROOT_OUT <- cfg('ROOT_OUT', file.path(getwd(), 'results'))
OUT_CSV   <- paste0(ROOT_OUT, "/bulkRNAseq_results")                 
OUT_FIG   <- paste0(ROOT_OUT, "/bulkRNAseq_figures") 

dir.create(OUT_CSV, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
EXPORT_DIR <- cfg("EXPORT_DIR", NULL) 

## 
## 1. load ----
## 
RDATA <- cfg('RDATA', file.path(getwd(), 'data', 'bulk_counts.RData'))
if (!exists("RAW_COUNTS2")) {
  if (!file.exists(RDATA)) stop("Input RData missing: ", RDATA,
                                "\nSet RDATA in config/config.R before running.")
  load(RDATA)
}
stopifnot(exists("RAW_COUNTS2"), exists("RAW_META2"))

RAW_META2 <- as.data.frame(RAW_META2)
stopifnot(all(c("sample", "type", "mouse") %in% colnames(RAW_META2)))

if (!identical(rownames(RAW_META2), as.character(RAW_META2$sample)))
  rownames(RAW_META2) <- as.character(RAW_META2$sample)

common <- intersect(rownames(RAW_META2), colnames(RAW_COUNTS2))
if (length(common) < ncol(RAW_COUNTS2))
  cat("NOTE:", ncol(RAW_COUNTS2) - length(common),
      "count columns have no metadata row and are dropped\n")
stopifnot(length(common) > 0)

RAW_META2   <- RAW_META2[common, , drop = FALSE]
RAW_COUNTS2 <- RAW_COUNTS2[, common, drop = FALSE]
stopifnot(identical(rownames(RAW_META2), colnames(RAW_COUNTS2)),
          !anyNA(RAW_META2$type), !anyNA(RAW_META2$mouse))

cat("=== loaded ===\n")
cat("genes:", nrow(RAW_COUNTS2), "| samples:", ncol(RAW_COUNTS2), "\n")
write.csv(RAW_META2, file.path(OUT_CSV, "input_metadata_aligned.csv"), row.names = TRUE)
print(table(region = RAW_META2$type))
print(table(mouse = RAW_META2$mouse, region = RAW_META2$type))

## 
## 2. two matrices: EC and Ctrl----
## 

EC_TYPES <- cfg("EC_TYPES", c("BEC", "CmEC", "ImEC", "IpEC"))
CONTROL_TYPES <- cfg("CONTROL_TYPES", c("Bt", "It"))

lv <- levels(factor(RAW_META2$type))
miss_ec  <- setdiff(EC_TYPES, lv)
miss_ctl <- setdiff(CONTROL_TYPES, lv)
if (length(miss_ec))
  stop("selection regions absent from metadata: ", paste(miss_ec, collapse = ", "),
       "\n  levels present: ", paste(lv, collapse = ", "))
if (length(miss_ctl))
  cat("NOTE: control regions absent:", paste(miss_ctl, collapse = ", "),
      "-- the figure will have no control columns\n")

sel <- rownames(RAW_META2)[as.character(RAW_META2$type) %in% EC_TYPES]
counts_EC <- RAW_COUNTS2[, sel, drop = FALSE]
meta_EC   <- RAW_META2[sel, , drop = FALSE]
meta_EC$type  <- factor(as.character(meta_EC$type), levels = EC_TYPES)
meta_EC$mouse <- factor(as.character(meta_EC$mouse))

all_t <- c(EC_TYPES, intersect(CONTROL_TYPES, lv))
alls  <- rownames(RAW_META2)[as.character(RAW_META2$type) %in% all_t]
counts_ALL <- RAW_COUNTS2[, alls, drop = FALSE]
meta_ALL   <- RAW_META2[alls, , drop = FALSE]
meta_ALL$type  <- factor(as.character(meta_ALL$type), levels = all_t)
meta_ALL$mouse <- factor(as.character(meta_ALL$mouse))

stopifnot(identical(rownames(meta_EC),  colnames(counts_EC)),
          identical(rownames(meta_ALL), colnames(counts_ALL)))

cat("\nselection matrix:", ncol(counts_EC), "samples /",
    nlevels(meta_EC$type), "regions | display matrix:",
    ncol(counts_ALL), "samples /", nlevels(meta_ALL$type), "regions\n")

## Counts must be integers for any NB model and for the voom
if (!all(counts_EC == round(counts_EC), na.rm = TRUE))
  cat("NOTE: non-integer counts detected (estimated counts?);",
      "rounding for the DE model\n")

##
## 3. knobs ----
##


LFC_MIN    <- cfg("LFC_MIN", 0.5)     # worst-case log2 FC floor 0.5
ALPHA      <- cfg("ALPHA", 0.05)    # BH-adjusted p-value threshold
N_TOP      <- cfg("N_TOP", 30)      # genes per region in the panel
MIN_DET    <- cfg("MIN_DET", 10)      # count floor for detection
NPERM      <- cfg("NPERM", 1000)    # permutation null for the accuracy number
DROP_MOUSE <- cfg("DROP_MOUSE", NULL)    # optional mouse ID to exclude
PAIRS_IN   <- cfg("PAIRS_IN", NULL)    # optional precomputed contrasts

## Merge regions the data cannot separate.
MERGE <- cfg("MERGE", list(mEC = c("CmEC", "ImEC")))

## 4. Figure parms ----
TOP_N        <- cfg("TOP_N", 5)                # genes per region shown
CONTROLS     <- cfg("CONTROLS", c("Bt", "It"))    # displayed, never used for selection
SHOW         <- cfg("SHOW", "samples")        # "samples" or "regions"
MARKER_SCALE <- cfg("MARKER_SCALE", "auto")           # "auto" | "z" | "absolute"
endothelial_markers <- cfg("endothelial_markers", c("Pecam1", "Cdh5", "Kdr", "Emcn", "Vwf"))


### design notes ----
##
## 1. Bt and It were NOT used to select genes. 
## 
## 2. z-scores are centred and scaled on the SELECTION REGIONS
##    ONLY, then the control samples are projected onto that
##    same scale. 
## 
## 3. The endothelial markers get their own scale and
##    colour map, never the signature z scale, because they
##    answer if a population endothelial.
##
##    Guard (MARKER_SCALE, default "auto"): z is used only when
##    depletion / within-EC scatter >= 2. Below that -- and
##    whenever no control region is present -- there is no real
##    contrast to scale against, z would amplify scatter between
##    flat EC columns, and absolute logCPM is the honest scale.
##    Override with MARKER_SCALE_IN <- "z" | "absolute".

## 
## 5. Selection, classification, diagnostics ----
## 
### 0. inputs ----
stopifnot(exists("counts_EC"), exists("meta_EC"))
cts <- counts_EC; meta <- meta_EC
stopifnot(all(colnames(cts) == rownames(meta)))
meta$type  <- droplevels(factor(meta$type))
meta$mouse <- droplevels(factor(meta$mouse))
GROUPS <- levels(meta$type)

if (!is.null(DROP_MOUSE)) {
  k <- as.character(meta$mouse) != DROP_MOUSE
  cts <- cts[, k]; meta <- droplevels(meta[k, ])
  cat("dropped mouse", DROP_MOUSE, "->", ncol(cts), "samples\n")
}
if (!is.null(MERGE)) {
  lv <- as.character(meta$type)
  for (nm in names(MERGE)) lv[lv %in% MERGE[[nm]]] <- nm
  meta$type <- droplevels(factor(lv))
  GROUPS <- levels(meta$type)
  cat("merged:", paste(names(MERGE), sapply(MERGE, paste, collapse = "+"),
                       sep = " <- ", collapse = "; "), "-> groups:",
      paste(GROUPS, collapse = ", "), "\n")
  cat("  NOTE: group sizes now unequal, so the relevant baseline is the\n",
      " MAJORITY class, not 1/k. Section 4 reports both.\n")
}

## 
### 1. all pairwise contrasts ----
## Every ordered pair (g, r). Blocking on mouse, since the data
## is paired across regions within animal.

pairwise <- function(cts, meta, groups) {
  keep <- rowSums(cts >= MIN_DET) >= min(table(meta$type))
  dge  <- calcNormFactors(DGEList(cts[keep, ]))
  has_mouse <- nlevels(meta$mouse) > 1 &&
    !is.null(tryCatch(model.matrix(~0 + type + mouse, meta), 
                      error = function(e) NULL))
  frm <- if (has_mouse) ~0 + type + mouse else ~0 + type
  des <- model.matrix(frm, meta)
  colnames(des) <- sub("^type", "", colnames(des))
  fit <- lmFit(voom(dge, des), des)
  gg  <- rownames(dge)
  L <- P <- array(NA_real_, c(length(gg), length(groups), 
                              length(groups)),
                  dimnames = list(gg, groups, groups))
  for (g in groups) for (r in setdiff(groups, g)) {
    cv <- setNames(rep(0, ncol(des)), colnames(des))
    cv[g] <- 1; cv[r] <- -1
    tt <- topTable(eBayes(contrasts.fit(fit, cv)), 
                   number = Inf, sort.by = "none")
    L[, g, r] <- tt$logFC; P[, g, r] <- tt$adj.P.Val
  }
  list(LFC = L, PADJ = P, genes = gg)
}

# PAIRS_IN is provided by config (NULL by default).
pw <- if (is.null(PAIRS_IN)) pairwise(cts, meta, GROUPS) else PAIRS_IN
LFC <- pw$LFC; PADJ <- pw$PADJ; genes <- pw$genes

cat("=== 1. pairwise contrasts ===\n")
cat("genes tested:", length(genes), "| contrasts:",
    length(GROUPS) * (length(GROUPS) - 1), "\n")
pc <- do.call(rbind, lapply(GROUPS, function(g)
  do.call(rbind, lapply(setdiff(GROUPS, g), function(r) data.frame(
    region = g, versus = r,
    n_sig = sum(PADJ[, g, r] < ALPHA & abs(LFC[, g, r]) >= LFC_MIN, na.rm = TRUE),
    n_up  = sum(PADJ[, g, r] < ALPHA & LFC[, g, r] >= LFC_MIN, na.rm = TRUE))))))
cat("\nsignificant genes per pairwise contrast (n_up = higher in `region`):\n")
print(pc[order(pc$region, pc$n_up), ], row.names = FALSE)

## 
### 2. worst-case separation score ----
## 
## min LFC and max padj across the other regions. Also the
## SECOND-worst, for the relaxed tier below.

wc <- lapply(GROUPS, function(g) {
  rv <- setdiff(GROUPS, g)
  lf <- LFC[, g, rv, drop = FALSE][, 1, ]; pj <- PADJ[, g, rv, drop = FALSE][, 1, ]
  ord <- t(apply(lf, 1, order))                       # rivals, hardest first
  hard <- rv[ord[, 1]]; second <- rv[ord[, 2]]
  data.frame(
    gene = genes, region = g,
    min_lfc     = apply(lf, 1, min, na.rm = TRUE),
    worst_padj  = apply(pj, 1, max, na.rm = TRUE),
    lfc2        = vapply(seq_len(nrow(lf)), function(i) lf[i, ord[i, 2]], numeric(1)),
    padj2       = vapply(seq_len(nrow(pj)), function(i)
      max(pj[i, ord[i, -1]], na.rm = TRUE), numeric(1)),
    hardest     = hard, second_hardest = second,
    row.names = NULL, stringsAsFactors = FALSE)
}) %>% bind_rows()

wc$tier <- with(wc, ifelse(min_lfc >= LFC_MIN & worst_padj < ALPHA, 1L,
                           ifelse(lfc2 >= LFC_MIN & padj2 < ALPHA, 2L, 
                                  NA_integer_)))

cat("\n=== 2. tier yield per region ===\n")
ty <- wc %>% filter(!is.na(tier)) %>% count(region, tier) %>%
  tidyr::pivot_wider(names_from = tier, values_from = n, values_fill = 0,
                     names_prefix = "tier")
print(as.data.frame(ty), row.names = FALSE)
cat("\n  tier 1 = separates from ALL other regions\n",
    " tier 2 = separates from all but its single closest rival\n")

## Which rival blocks each region most often -- the number that
## tells you whether two of the regions are simply not distinct 
## in their expression profile.
cat("\nblocking rival for tier-2 genes:\n")
print(as.data.frame(wc %>% filter(tier == 2L) %>% count(region, hardest) %>%
                      group_by(region) %>% slice_max(n, n = 2) %>% ungroup()), 
      row.names = FALSE)

## 
### 3. assemble the panel ----
## 
## Tier 1 first, ranked by worst-case margin. Top up from tier 2
## only if a region is short, and label it so the distinction
## survives into the figure and the paper.

## A gene can qualify as tier 2 for BOTH members of a
## non-distinct pair -- each region's only blocker is the other.
## Those are shared-module genes: the exact overlap this
## selection exists to prevent. Drop them rather than assign
## them twice, and count them, because a large count IS the
## finding that two regions are not separable.
t2_all <- wc %>% filter(tier == 2L)
recip  <- t2_all$gene[duplicated(t2_all$gene)] %>% unique()
if (length(recip))
  cat("\n  NOTE:", length(recip), 
      "genes qualify as tier 2 for >1 region",
      "(shared modules) -- excluded\n")
wc$tier[wc$tier == 2L & wc$gene %in% recip] <- NA_integer_

panel <- lapply(GROUPS, function(g) {
  s  <- wc %>% filter(region == g, !is.na(tier))
  t1 <- s %>% filter(tier == 1L) %>% arrange(desc(min_lfc)) %>% head(N_TOP)
  need <- N_TOP - nrow(t1)
  t2 <- if (need > 0) s %>% filter(tier == 2L) %>% arrange(desc(lfc2)) %>% head(need)
  else s[0, ]
  bind_rows(t1, t2)
}) %>% bind_rows() %>%
  mutate(score = ifelse(tier == 1L, min_lfc, lfc2)) %>%
  arrange(region, tier, desc(score))

cat("\n=== 3. panel composition ===\n")
if (!nrow(panel)) {
  cat("  NO GENES PASS. Nothing separates any region at",
      "LFC >=", LFC_MIN, "and FDR <", ALPHA, ".\n",
      " Lower LFC_MIN, or accept that these regions are not\n",
      " distinguishable at this depth. Section 4 still runs and\n",
      " will report near-chance accuracy, which is the answer.\n")
} else {
  print(as.data.frame(panel %>% count(region, tier)), row.names = FALSE)
  stopifnot(!any(duplicated(panel$gene)))   # non-overlap is structural
  cat("\nduplicate genes across regions:", sum(duplicated(panel$gene)),
      "(zero by construction)\n")
}

if (nrow(panel)) {
  cat("\ntop 5 per region (score = worst-case log2 FC):\n")
  print(as.data.frame(panel %>% group_by(region) %>% slice_head(n = 5) %>%
                        ungroup() %>% dplyr::select(region, 
                                                    gene, 
                                                    tier, 
                                                    score, 
                                                    hardest) %>%
                        mutate(score = round(score, 2))), row.names = FALSE)
}

#### outputs----
write.csv(panel, file.path(OUT_CSV, "region_panel.csv"), row.names = FALSE)
write.csv(wc, file.path(OUT_CSV, "worst_case_scores.csv"), row.names = FALSE)

##
### 4. Evaluation (LOMO)----
##
## LOLO: Leave-one-MOUSE-out. Gene selection is
## recalculated inside every fold otherwise the held-out sample
## helped choose its own classifier and the accuracy is inflated.

lcpm_of <- function(cc) cpm(calcNormFactors(DGEList(cc)), 
                            log = TRUE, prior.count = 2)

classify <- function(train_cts, train_meta, test_cts, n_per = N_TOP) {
  p  <- pairwise(train_cts, train_meta, GROUPS)
  gs <- intersect(p$genes, rownames(test_cts))
  sel   <- unlist(lapply(GROUPS, function(g) {
      rv  <- setdiff(GROUPS, g)
      lf  <- p$LFC[gs, g, rv, drop = FALSE][, 1, ]
      pj  <- p$PADJ[gs, g, rv, drop = FALSE][, 1, ]
      m   <- apply(lf, 1, min, na.rm = TRUE) 
      q   <- apply(pj, 1, max, na.rm = TRUE)
      ok <- which(q < ALPHA & m >= LFC_MIN)
    ## FALLBACK: nothing passed. Take the best-ranked genes anyway so
    ## the fold still runs -- but record it. A region whose panel is
    ## always forced has no signature, and its accuracy is not
    ## evidence about that region.
    if (!length(ok)) {
      FORCED[[g]] <<- c(FORCED[[g]], 1L)
      ok <- order(m, decreasing = TRUE)[seq_len(min(n_per, length(m)))]
    }
    names(sort(m[ok], decreasing = TRUE))[seq_len(min(n_per, length(ok)))]
  }))
  sel   <- unique(na.omit(sel))
  lall  <- lcpm_of(cbind(train_cts[gs, ], test_cts[gs, , drop = FALSE]))
  mu    <- rowMeans(lall[sel, , drop = FALSE]) 
  sd_   <- apply(lall[sel, , drop = FALSE], 1, sd)
  sd_[sd_ == 0] <- 1
  z     <- (lall[sel, , drop = FALSE] - mu) / sd_
  ztr   <- z[, colnames(train_cts), drop = FALSE]
  zte   <- z[, colnames(test_cts), drop = FALSE]
  cen <- vapply(GROUPS, function(g)
    rowMeans(ztr[, train_meta$type == g, drop = FALSE]), numeric(nrow(ztr)))
  apply(zte, 2, function(x) GROUPS[which.max(cor(x, cen))])
}

mice <- levels(meta$mouse)
cat("\n=== 4. leave-one-mouse-out classification ===\n")
cat("folds:", length(mice), "| selection redone inside each fold\n")
FORCED <- setNames(vector("list", length(GROUPS)), GROUPS)
pred <- character(0); truth <- character(0); sname <- character(0)
pmouse <- character(0)
for (m in mice) {
  te <- as.character(meta$mouse) == m
  if (sum(!te) < 2 * length(GROUPS)) { 
    cat("  skip", m, "(too few left)\n"); next
    }
  tm <- droplevels(meta[!te, ]); if (nlevels(tm$type) < length(GROUPS)) next
  p <- classify(cts[, !te], tm, cts[, te, drop = FALSE])
  pred <- c(pred, p); truth <- c(truth, as.character(meta$type[te]))
  sname <- c(sname, rownames(meta)[te]); pmouse <- c(pmouse, rep(m, sum(te)))
}
cm <- table(true = factor(truth, GROUPS), 
            predicted = factor(pred, GROUPS))

cat("\nconfusion matrix:\n"); print(cm)
acc <- data.frame(region = GROUPS,
                  n = as.integer(rowSums(cm)),
                  correct = as.integer(diag(cm)),
                  accuracy = round(diag(cm) / pmax(rowSums(cm), 1), 2))

acc$most_confused_with <- vapply(GROUPS, function(g) {
  o <- cm[g, setdiff(GROUPS, g)]
  if (all(o == 0)) "-" else names(which.max(o)) }, character(1))
acc$panel_forced <- vapply(GROUPS, function(g) length(FORCED[[g]]), integer(1))
acc$never_predicted <- colSums(cm)[GROUPS] == 0
cat("\nper-region accuracy:\n")
print(acc, row.names = FALSE)

## Baselines. With unequal group sizes 1/k is the WRONG baseline --
## always predicting the largest class beats it. Report both, plus a
## permutation null that respects the mouse pairing (labels shuffled
## WITHIN mouse), which is the only one that accounts for the fact
## that n is small and the folds are not independent.
nk    <- rowSums(cm)
p_maj <- max(nk) / sum(nk)
obs_a <- sum(diag(cm)) / sum(cm)
set.seed(1)
ix_by_mouse <- split(seq_along(truth), pmouse)
perm <- replicate(NPERM, {
  hit <- 0L
  for (ix in ix_by_mouse) hit <- hit + sum(sample(pred[ix]) == truth[ix])
  hit / length(truth)
})
cat("\noverall:", round(obs_a, 3),
    "| uniform 1/k:", round(1 / length(GROUPS), 3),
    "| majority class:", round(p_maj, 3), "\n")
cat("permutation null (labels shuffled within mouse,", NPERM, "draws):",
    "mean", round(mean(perm), 3), "| p =",
    signif((1 + sum(perm >= obs_a)) / (NPERM + 1), 3), "\n")

if (any(acc$never_predicted))
  cat("\n  WARNING:", paste(GROUPS[acc$never_predicted], collapse = ", "),
      "was never predicted in any fold.\n",
      " Its accuracy of 0 is not 'poorly identified' -- the classifier\n",
      " has no usable panel for it, so samples are absorbed by its\n",
      " nearest neighbour. Check `panel_forced`: if that region's panel\n",
      " was forced in every fold, the region has no signature at this\n",
      " depth and MERGE is the honest move.\n")
if (any(acc$panel_forced > 0))
  cat("\n  panel forced (no gene passed; best-ranked used instead):",
      paste(sprintf("%s %d/%d folds", GROUPS[acc$panel_forced > 0],
                    acc$panel_forced[acc$panel_forced > 0], length(mice)),
            collapse = "; "), "\n")

#### outputs----
write.csv(data.frame(sample = sname, true = truth, predicted = pred),
          file.path(OUT_CSV, "loo_predictions.csv"), row.names = FALSE)

write.csv(acc, file.path(OUT_CSV, "classification_accuracy.csv"), 
          row.names = FALSE)

## 
## 6. Figures ----
## 
lk <- lcpm_of(cts)
gm <- vapply(GROUPS, function(g)
  rowMeans(lk[, meta$type == g, drop = FALSE]), numeric(nrow(lk)))
pg <- intersect(panel$gene, rownames(gm))

if (length(pg) >= 2) {
  z  <- t(scale(t(gm[pg, , drop = FALSE])))
  
  hd <- as.data.frame(as.table(z)); names(hd) <- c("gene", "region", "z")
  hd$region <- factor(hd$region, levels = GROUPS)
  ord <- panel %>% filter(gene %in% pg) %>% arrange(region, desc(score))
  hd$gene <- factor(hd$gene, levels = rev(ord$gene))
  hd$block <- panel$region[match(as.character(hd$gene), panel$gene)]
  hd$tier <- panel$tier[match(as.character(hd$gene), panel$gene)]
  tier_strip <- hd %>% distinct(gene, tier) %>% mutate(x = "tier")
  n_t2 <- sum(panel$tier == 2L)
  
  p_heat <- ggplot(hd, aes(region, gene, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2C6E8E", mid = "grey96", high = "#B23A34",
                         midpoint = 0, name = "z") +
    ggnewscale::new_scale_fill() +
    geom_tile(data = tier_strip, aes(x = x, y = gene, fill = factor(tier)),
              inherit.aes = FALSE, width = .5) +
    scale_fill_manual(values = c(`1` = "grey30", `2` = "#E8A33D"),
                      name = "tier", 
                      labels = c("1: beats all", "2: one exception")) +
    labs(title = "Region panel, worst-case selection",
         subtitle = paste0(nrow(ord), " genes: ", nrow(ord) - n_t2,
                           " beat every other region by >= 2^", LFC_MIN,
                           ";\n", n_t2, " (orange) beat all but one named rival"),
         x = NULL, y = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_text(size = 4.5),
          axis.ticks.y = element_blank(),
          legend.key.width = unit(7, "pt"),
          legend.key.height = unit(7, "pt"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6.5),
          plot.title = element_text(size = 9.5),
          plot.subtitle = element_text(size = 7, colour = "grey35"))
  ggsave(file.path(OUT_FIG, "region_panel_heatmap.png"), p_heat,
         width = 3.6, height = 7.2, dpi = 300, limitsize = FALSE)
} else cat("\n  (no heatmap: panel is empty)\n")

cmd <- as.data.frame(as.table(cm))
p_cm <- ggplot(cmd, aes(predicted, true, fill = Freq)) +
  geom_tile(colour = "white", linewidth = .6) +
  geom_text(aes(label = Freq), size = 3) +
  scale_fill_gradient(low = "grey96", high = "#2C6E8E", guide = "none") +
  labs(title = "Leave-one-mouse-out confusion",
       subtitle = paste0("overall ", round(100 * sum(diag(cm)) / sum(cm)),
                         "% vs ", round(100 / length(GROUPS)), "% chance"),
       x = "predicted", y = "true") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9.5),
        plot.subtitle = element_text(size = 7.6, colour = "grey35"))
ggsave(file.path(OUT_FIG, "confusion_matrix.png"), p_cm,
       width = 3.8, height = 3.2, dpi = 300)

saveRDS(list(panel = panel, wc = wc, cm = cm, acc = acc,
             LFC = LFC, PADJ = PADJ),
        file.path(OUT_CSV, "region_classifier.rds"))

## 
## 5b. export shapes
## 
## Long `panel` is the source of truth. These are convenience
## shapes -- but note that both of them drop `tier`, and a
## tier-2 gene is NOT interchangeable with a tier-1 gene: it is
## equally high in the region named in `hardest`. Handing a
## flattened gene list to a scoring function reintroduces
## exactly the overlap this script exists to remove. So the
## wide export carries a tier mark, and the list is also
## provided tier-1-only.

if (nrow(panel)) {
  ## named list, for scoring / GSEA / loops
  sig_list <- split(panel$gene, panel$region)
  t1 <- panel[panel$tier == 1L, ]
  ## split on a factor with only the OBSERVED tier-1 regions as levels, so
  ## a region with no tier-1 gene is absent from the list rather than
  ## present-and-empty -- silently scoring an empty signature is worse
  ## than getting a "no such element" error.
  sig_list_t1 <- if (nrow(t1))
    split(t1$gene, factor(t1$region, levels = intersect(GROUPS, unique(t1$region))))
  else list()
  
  ## padded wide table, for a supplementary file. pivot_wider needs an
  ## explicit within-region index -- without id_cols every row is its
  ## own group and you get one row of list-columns, which unlist()
  ## then flattens back to the vector you started with.
  wide <- panel %>%
    mutate(label = ifelse(tier == 2L, paste0(gene, " (=", hardest, ")"), gene)) %>%
    group_by(region) %>% mutate(rank = row_number()) %>% ungroup() %>%
    tidyr::pivot_wider(id_cols = rank, names_from = region,
                       values_from = label) %>%
    arrange(rank) %>% dplyr::select(-rank) %>%
    mutate(across(everything(), ~ifelse(is.na(.), "", .)))
  
  write.csv(wide, file.path(OUT_CSV, "signature_wide.csv"), row.names = FALSE)
  saveRDS(list(long = panel, list_all = sig_list, list_tier1 = sig_list_t1),
          file.path(OUT_CSV, "signature_shapes.rds"))
  
  cat("\n=== 5b. exports ===\n")
  cat("signature_wide.csv:", nrow(wide), "x", ncol(wide),
      "- tier-2 genes marked '(=rival)'\n")
  cat("sig_list:   ", paste(sprintf("%s=%d", names(sig_list),
                                    lengths(sig_list)), collapse = " "), "\n")
  cat("sig_list_t1:", paste(sprintf("%s=%d", names(sig_list_t1),
                                    lengths(sig_list_t1)), collapse = " "), "<- use for scoring\n")
  gone <- setdiff(GROUPS, names(sig_list_t1))
  if (length(gone))
    cat("  NOTE:", paste(gone, collapse = ", "),
        "has NO tier-1 gene, so it is ABSENT from\n",
        " sig_list_t1 rather than present-and-empty. A tier-2-only\n",
        " region has no standalone signature: every gene it has is\n",
        " shared with the rival named in panel$hardest.\n")
}

## 
## Extra. how to read this----
## 
cat("\n=== 6. how to read it ===\n")
cat("* tier 1 genes are the ones to report. Non-overlap is\n",
    "  structural, so no de-duplication and no tie-breaking.\n",
    "* tier 2 genes separate a region from everything except one\n",
    "  rival, named in the `hardest` column. They are still\n",
    "  useful, but state the exception.\n",
    "* section 4 accuracy is the honest answer to \"does this\n",
    "  identify the region\". A region at 100% with only 3 tier-1\n",
    "  genes is identified fine -- gene COUNT is not the goal.\n",
    "* a region that classifies poorly and whose tier-2 blocker\n",
    "  is always the same partner is not a separate region at\n",
    "  this depth. Merging them is a legitimate result.\n")

## 
## 7. final signature heatmap ----
## 
## ---- inputs ----
## counts_ALL / meta_ALL: ALL SIX regions (unlike the classifier,
## which was run on the four). meta_ALL needs $type and $mouse.
stopifnot(exists("counts_ALL"), exists("meta_ALL"), exists("panel"))

cts <- counts_ALL; meta <- as.data.frame(meta_ALL)
stopifnot(all(colnames(cts) == rownames(meta)))
meta$type  <- factor(as.character(meta$type))
meta$mouse <- factor(as.character(meta$mouse))

## Re-apply MERGE here. This section rebuilds labels from meta_ALL,
## which is the RAW six-region metadata -- it never went through the
## merge in section 3. Without this the panel says "mEC" while SEL
## says "CmEC","ImEC", so `filter(region %in% SEL)` below silently
## drops every merged gene and the figure loses whole row blocks.
if (!is.null(MERGE)) {
  lv <- as.character(meta$type)
  for (nm in names(MERGE)) lv[lv %in% MERGE[[nm]]] <- nm
  meta$type <- factor(lv)
  cat("figure: merged", paste(names(MERGE),
                              sapply(MERGE, paste, collapse = "+"), sep = " <- ", collapse = "; "),
      "\n")
}

SEL <- setdiff(levels(meta$type), CONTROLS)     # selection regions
CTL <- intersect(CONTROLS, levels(meta$type))   # controls actually present
if (!length(CTL))
  cat("NOTE: none of", paste(CONTROLS, collapse = "/"), "found in meta_ALL$type.",
      "Levels are:", paste(levels(meta$type), collapse = ", "), "\n")
meta$type <- factor(as.character(meta$type), levels = c(SEL, CTL))
meta$role <- ifelse(as.character(meta$type) %in% CTL,
                    "control (not used for selection)", "selection regions")
meta$role <- factor(meta$role,
                    levels = c("selection regions", "control (not used for selection)"))

cat("=== inputs ===\n")
cat("selection regions:", paste(SEL, collapse = ", "),
    "| controls:", ifelse(length(CTL), paste(CTL, collapse = ", "), "none"), "\n")
print(table(region = meta$type))

## 
## 1. gene blocks ----
## 
## Top TOP_N per region, tier 1 first. A region whose top-5 is
## partly tier 2 is flagged in the row label, because a tier-2
## gene is equally high in the rival named in `hardest`.

.orphan <- setdiff(unique(as.character(panel$region)), SEL)
if (length(.orphan))
  stop("panel has regions absent from the figure's SEL: ",
       paste(.orphan, collapse = ", "),
       "\n  SEL is: ", paste(SEL, collapse = ", "),
       "\n  usually means MERGE was applied to the selection stage but",
       " not to meta_ALL.")

top <- panel %>%
  filter(region %in% SEL) %>%
  arrange(region, tier, desc(score)) %>%
  group_by(region) %>% 
  slice_head(n = TOP_N) %>% 
  ungroup() %>%
  mutate(block = as.character(region),
         lab = ifelse(tier == 2L, 
                      paste0(gene, " *"), gene))

## Count against ALL selection regions, not just those present in `top`.
## A region with ZERO genes is absent from count() entirely, so a naive
## filter(n < TOP_N) silently misses the worst case -- and the figure
## would then show that region's COLUMN with no row block, which a
## reader can easily fail to notice.
navail <- setNames(rep(0L, length(SEL)), SEL)
tb <- table(as.character(top$region)); navail[names(tb)] <- as.integer(tb)
short <- navail[navail < TOP_N]
EMPTY <- names(navail)[navail == 0L]

if (length(short))
  cat("\nNOTE: fewer than", TOP_N, "genes available for:",
      paste(sprintf("%s (%d)", names(short), short), collapse = ", "), "\n")
if (length(EMPTY))
  cat("*** ", paste(EMPTY, collapse = ", "),
      "contributed NO genes to the panel. The column is still\n",
      "    drawn (its samples exist) but there is no row block for it.",
      "That is\n    the finding: no gene distinguishes it from every other region.\n")

miss <- setdiff(endothelial_markers, rownames(cts))
if (length(miss))
  cat("NOTE: markers absent from the count matrix:",
      paste(miss, collapse = ", "), "\n")
em <- intersect(endothelial_markers, rownames(cts))

n_t2 <- sum(top$tier == 2L)
cat("\npanel for figure:", nrow(top), "signature genes +", length(em),
    "markers |", n_t2, "are tier 2 (marked *)\n")

## 
## 2. expression matrices ----
## 
keep <- rowSums(cts >= MIN_DET) >= min(table(meta$type))
keep[rownames(cts) %in% c(top$gene, em)] <- TRUE
lcpm <- cpm(calcNormFactors(DGEList(cts[keep, ])), 
            log = TRUE, prior.count = 2)

sel_cols <- rownames(meta)[!as.character(meta$type) %in% CTL]

## z on selection samples only, controls projected onto it 
zscore_proj <- function(m, ref_cols) {
  mu <- rowMeans(m[, ref_cols, drop = FALSE])
  sd <- apply(m[, ref_cols, drop = FALSE], 1, sd) 
  sd[sd == 0 | is.na(sd)] <- 1
  (m - mu) / sd
}

if (SHOW == "regions") {
  agg <- vapply(levels(meta$type), function(g)
    rowMeans(lcpm[, meta$type == g, drop = FALSE]), numeric(nrow(lcpm)))
  mat  <- agg
  cmeta <- data.frame(col = colnames(agg), 
                      type = factor(colnames(agg),
                                    levels = levels(meta$type)), 
                      role = ifelse(colnames(agg) %in% CTL,
                                                                                                    "control (not used for selection)", "selection regions"))
  ref <- setdiff(colnames(agg), CTL)
} else {
  ord <- order(meta$type, meta$mouse)
  mat <- lcpm[, ord, drop = FALSE]
  cmeta <- data.frame(col = rownames(meta)[ord], type = meta$type[ord],
                      role = meta$role[ord])
  ref <- intersect(sel_cols, colnames(mat))
}
cmeta$col <- factor(cmeta$col, levels = cmeta$col)
cmeta$role <- factor(as.character(cmeta$role),
                     levels = c("selection regions", "control (not used for selection)"))

zm <- zscore_proj(mat, ref)

long_sig <- as.data.frame(as.table(zm[top$gene, , drop = FALSE])) %>%
  setNames(c("gene", "col", "z")) %>%
  left_join(top %>% dplyr::select(gene, block, lab, tier), by = "gene") %>%
  left_join(cmeta, by = "col")

## ---- marker scale (see header decision 3) ----------------
## The markers are NOT expected to be flat across all columns: they are
## flat across the EC regions and DEPLETED in the tumour-only controls.
## That depletion is the positive control, so the scale has to show it.
##
## z must therefore be computed across ALL columns, never selection-only:
## within the EC regions the markers are flat, so a selection-only SD is a
## noise SD, and dividing the real depletion by it sends the control tiles
## to z -7.9 .. -15.1 (measured) -- off any usable colour scale. Including
## the controls puts the genuine EC-vs-tumour variance in the denominator.
## MARKER_SCALE is set in section 3.
ctype <- setNames(as.character(cmeta$type), cmeta$col)
ec_c  <- colnames(mat)[!ctype[colnames(mat)] %in% CTL]
ct_c  <- colnames(mat)[ ctype[colnames(mat)] %in% CTL]

em_z <- NULL; mk_guard <- NULL
if (length(em)) {
  mm  <- mat[em, , drop = FALSE]
  sdm <- apply(mm, 1, sd)
  em_z <- sweep(mm, 1, rowMeans(mm), "-") / ifelse(sdm > 0, sdm, NA_real_)
  if (length(ct_c)) {
    gap    <- mean(rowMeans(em_z[, ec_c, drop = FALSE]) -
                     rowMeans(em_z[, ct_c, drop = FALSE]))
    spread <- mean(apply(em_z[, ec_c, drop = FALSE], 1, sd))
    mk_guard <- c(gap = gap, spread = spread, ratio = gap / spread)
  }
}
## Guard: only use z if the depletion is large relative to within-EC
## scatter. If it is not, z is amplifying noise across the EC columns and
## absolute logCPM is the honest scale -- which is also the only option
## when no control is present, since then there is no real variance to
## scale against.
USE_Z_MK <- length(em) > 0 &&
  (MARKER_SCALE == "z" ||
     (MARKER_SCALE == "auto" && !is.null(mk_guard) && mk_guard["ratio"] >= 2))
MK_NAME <- if (USE_Z_MK) "z\n(all regions)" else "logCPM\n(absolute)"

long_em <- as.data.frame(as.table(
  if (USE_Z_MK) em_z else mat[em, , drop = FALSE])) %>%
  setNames(c("gene", "col", "mkval")) %>%
  mutate(block = "endothelial\nmarkers", lab = as.character(gene)) %>%
  left_join(cmeta, by = "col")

## row order: signature blocks in region order, then the QC block
row_lv <- c(rev(top$lab), rev(em))
long_sig$lab <- factor(long_sig$lab, levels = row_lv)
long_em$lab  <- factor(long_em$lab,  levels = row_lv)
blk_lv <- c(SEL, "endothelial\nmarkers")
long_sig$block <- factor(long_sig$block, levels = blk_lv)
long_em$block  <- factor(long_em$block,  levels = blk_lv)

## report the control projection: how far outside the panel range?
if (length(CTL)) {
  zc <- long_sig %>% filter(as.character(type) %in% CTL)
  cat("\ncontrol columns, z relative to the selection regions:",
      "median", round(median(zc$z), 2),
      "| range", paste(round(range(zc$z), 2), collapse = " to "), "\n")
}

## The QC claim, stated numerically rather than left to the eye. The
## expectation is directional: flat and high across the EC regions, LOW in
## the tumour-only controls. Both halves are testable, so test both and
## name any marker that fails, rather than printing one median per region
## and leaving the reader to eyeball a contrast.
if (length(em)) {
  emq <- long_em %>% group_by(type) %>%
    summarise(median = round(median(mkval), 2), .groups = "drop")
  cat("\nendothelial markers, median", ifelse(USE_Z_MK, "z", "logCPM"),
      "by region:\n")
  print(as.data.frame(emq), row.names = FALSE)
  
  if (length(ct_c)) {
    lfc <- rowMeans(mat[em, ec_c, drop = FALSE]) -
      rowMeans(mat[em, ct_c, drop = FALSE])
    ecsd <- apply(mat[em, ec_c, drop = FALSE], 1, sd)
    mkt <- data.frame(marker = em,
                      EC_mean  = round(rowMeans(mat[em, ec_c, drop = FALSE]), 2),
                      ctl_mean = round(rowMeans(mat[em, ct_c, drop = FALSE]), 2),
                      logFC_EC_vs_ctl = round(lfc, 2),
                      within_EC_sd = round(ecsd, 2))
    cat("\nEC regions vs tumour-only controls (", paste(CTL, collapse = "/"),
        "), logCPM:\n", sep = "")
    print(mkt, row.names = FALSE)
    weak <- em[lfc < 1]
    if (length(weak))
      cat("*** weak/absent depletion (logFC < 1):", paste(weak, collapse = ", "),
          "\n    these do not discriminate EC from tumour here, so they carry",
          "no QC\n    weight -- consider dropping them from the block.\n")
    if (!is.null(mk_guard))
      cat("depletion / within-EC scatter =", round(mk_guard["ratio"], 1),
          "->", ifelse(USE_Z_MK, "z scale (real contrast to scale against)",
                       "absolute scale (z would amplify EC-column noise)"), "\n")
  } else {
    cat("  no tumour-only control present, so the EC-vs-tumour contrast",
        "cannot be\n  shown; marker block is on absolute logCPM.\n")
  }
}

## 
## 3. the figure ----
## 
zlim <- c(-2.2, 2.2)
long_sig$z_cl <- pmax(pmin(long_sig$z, zlim[2]), zlim[1])
n_out <- sum(long_sig$z > zlim[2] | long_sig$z < zlim[1])

if (n_out > 0)
  cat("\nNOTE:", n_out, "of", nrow(long_sig), "signature tiles",
      sprintf("(%.1f%%)", 100 * n_out / nrow(long_sig)),
      "exceed z =", zlim[2], "and are drawn clipped.\n",
      "  observed z range:", paste(round(range(long_sig$z), 2), collapse = " to "),
      "\n")

base_theme <- theme_minimal(base_size = 8.4) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = .5, 
                                   size = 5.6),
        axis.text.y = element_text(size = 6.2),
        strip.text.y.left = element_text(angle = 0, 
                                         size = 6.6, 
                                         face = "bold"),
        strip.text.x = element_text(size = 6.8),
        strip.placement = "outside",
        panel.spacing.x = unit(.14, "lines"),
        panel.spacing.y = unit(.30, "lines"),
        legend.key.width = unit(.30, "cm"),
        legend.key.height = unit(.42, "cm"),
        legend.title = element_text(size = 6.8),
        legend.text = element_text(size = 6.2),
        plot.title = element_text(size = 9.6),
        plot.subtitle = element_text(size = 7.0, 
                                     colour = "grey35"),
        plot.caption = element_text(size = 6.1, 
                                    colour = "grey40", 
                                    hjust = 0))

p <- ggplot(long_sig, aes(col, lab, fill = z_cl)) +
  geom_tile() +
  scale_fill_gradient2(high = "darkblue", 
                       mid = "grey96", 
                       low = "#C8A21C",
                       midpoint = 0, 
                       limits = zlim, 
                       oob = scales::squish,
                       name = "z\n(vs selection\nregions)") +
  new_scale_fill() +
  geom_tile(data = long_em, aes(col, lab, fill = mkval), inherit.aes = FALSE) +
  (if (USE_Z_MK)
    scale_fill_gradient2(high = "darkgreen", 
                         mid = "grey96", 
                         low = "grey25",
                         midpoint = 0, 
                         name = MK_NAME)
   else
     scale_fill_viridis_c(option = "mako", 
                          direction = -1, 
                          name = MK_NAME)) +
  facet_grid(block ~ type, 
             scales = "free", 
             space = "free", 
             switch = "y") +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(axis.text.x = element_blank())


strip <- ggplot(cmeta, aes(col, 1, fill = role)) +
  geom_tile() +
  facet_grid(. ~ type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("selection regions" = "black",
                               "control (not used for selection)" = "grey70"), 
                    name = NULL) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8.4) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = .5, 
                                   size = 5.6),
        axis.text.y = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing.x = unit(.14, "lines"),
        legend.position = "bottom",
        legend.key.size = unit(.30, "cm"),
        legend.text = element_text(size = 6.2),
        legend.margin = margin(t = -2, b = 0))

pf <- (p / strip) +
  plot_layout(heights = c(nrow(top) + length(em), 1.05))

print(pf)

W <- if (SHOW == "regions") 5.6 else 2.5 + .105 * ncol(mat)
H <- 2.4 + .128 * (nrow(top) + length(em))

#### save figure signature ----
ggsave(file.path(OUT_FIG, "final_signature_heatmap.png"), pf,
       width = W, height = H, dpi = 300, limitsize = FALSE)
cat("\nsaved final_signature_heatmap.png (", round(W,1), "x", round(H,1), "in )\n")


if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file.path(OUT_FIG, "final_signature_heatmap.svg"),
                   width = W, height = H)
  print(pf); dev.off()
  cat("saved final_signature_heatmap.svg (same", round(W,1), "x", round(H,1), "in )\n")
} else {
  cat("NOTE: svglite not installed -- SVG skipped, PNG written.\n",
      "  install.packages(\"svglite\") if you want the vector version.\n")
}

#### save signature used in figure----
write.csv(top %>% dplyr::select(region, gene, tier, score, hardest),
          file.path(OUT_CSV, "figure_genes.csv"), row.names = FALSE)


## 
## 5. export csv ----
## 

wide <- panel %>%
  group_by(region) %>% mutate(rank = row_number()) %>% ungroup() %>%
  pivot_wider(id_cols = rank, names_from = region, values_from = gene) %>%
  arrange(rank) %>% dplyr::select(-rank) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

write.csv(wide, file.path(OUT_CSV, "ec_signatures_20_wide.csv"), row.names = FALSE)

if (!is.null(EXPORT_DIR)) {
  if (dir.exists(EXPORT_DIR)) {
    dir.create(EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)
    for (f in c("ec_signatures_20_wide.csv", "region_panel.csv",
                "figure_genes.csv", "classification_accuracy.csv")) {
      from <- file.path(OUT_CSV, f)
      if (file.exists(from)) file.copy(from, file.path(EXPORT_DIR, f), overwrite = TRUE)
    }
    from <- file.path(OUT_FIG, "final_signature_heatmap.png")
    if (file.exists(from))
      file.copy(from, file.path(EXPORT_DIR, basename(from)), overwrite = TRUE)
    cat("\ncopied deliverables to", EXPORT_DIR, "\n")
  } else cat("\nEXPORT_DIR not reachable, skipped:", EXPORT_DIR, "\n")
}

cat("\n=== done ===\noutputs in", normalizePath(OUT_CSV), "\n")

## 
## 8. GO enrichment and GSEA
## 
## Two different questions, deliberately both:
##
##   GSEA  -- uses EVERY tested gene, ranked. Asks "is this GO
##            term shifted up in region g relative to the other
##            regions?" No threshold, so it is not sensitive to
##            LFC_MIN / N_TOP. This is the one that reports NES.
##   ORA   -- hypergeometric test on the SELECTED panel genes
##            only. Asks "is this term over-represented among the
##            genes I put in the signature?" It inherits every
##            selection threshold, so read it as a description of
##            the panel, not as independent evidence.
##
## GSEA is primary. ORA is reported because it is what the panel
## itself is made of, and disagreement between the two is
## informative: a term that is ORA-significant but not GSEA-
## significant is riding on a handful of top genes.

## ---- knobs (the only place these are set) ----
RUN_GO        <- cfg("RUN_GO", TRUE)
GO_ORGDB      <- cfg("GO_ORGDB", "org.Mm.eg.db")   # mouse. Human: "org.Hs.eg.db"
GO_ONT        <- cfg("GO_ONT", "BP")             # "BP" | "MF" | "CC" or a vector
GO_MINSIZE    <- cfg("GO_MINSIZE", 15)               # AFTER restricting to the universe
GO_MAXSIZE    <- cfg("GO_MAXSIZE", 500)              # huge terms are uninformative
GO_RANK       <- cfg("GO_RANK", "mean_lfc")       # "mean_lfc" | "min_lfc"; see note
GO_NPERM      <- cfg("GO_NPERM", 10000)            # fgsea preliminary permutations
GO_TOP_PLOT   <- cfg("GO_TOP_PLOT", 5)                # terms per region in the figures
GO_JACCARD    <- cfg("GO_JACCARD", 0.5)              # leading-edge overlap -> redundant
GO_PADJ       <- cfg("GO_PADJ", 0.05)

## GO_RANK note. The selection statistic is `min_lfc`, the WORST
## case across rivals. That is right for picking genes but wrong
## as a GSEA ranking: a min over k-1 rivals is systematically
## shifted negative, so the down tail stops meaning "lower in
## this region" and NES for depleted terms becomes hard to read.
## `mean_lfc` -- the one-vs-rest average -- is signed symmetric
## and is the default. The script prints both distributions so
## you can see the shift rather than take my word for it.

if (!RUN_GO) {
  cat("\n=== 8. GO: skipped (RUN_GO = FALSE) ===\n")
} else {
  
  .go_ok <- all(vapply(c("fgsea", "AnnotationDbi", "GO.db", GO_ORGDB),
                       function(p) requireNamespace(p, quietly = TRUE), logical(1)))
  if (!.go_ok) {
    cat("\n=== 8. GO: SKIPPED -- missing packages ===\n")
    cat("  install with:\n",
        '  BiocManager::install(c("fgsea","org.Mm.eg.db","GO.db","AnnotationDbi"))\n')
  } else {
    
    suppressPackageStartupMessages({
      library(fgsea); library(AnnotationDbi); library(GO.db); library(GO_ORGDB, character.only = TRUE)
    })
    .odb <- get(GO_ORGDB)
    
    cat("\n=========================================================\n")
    cat("=== 8. GO enrichment and GSEA ===\n")
    cat("=========================================================\n")
    
    ## 
    ## 8.1 symbol -> Entrez ----
    ## 
    ## Two passes: official SYMBOL, then ALIAS for whatever failed.

    .map1 <- suppressMessages(AnnotationDbi::mapIds(
      .odb, keys = genes, keytype = "SYMBOL", column = "ENTREZID",
      multiVals = "first"))
    .todo <- genes[is.na(.map1)]
    if (length(.todo)) {
      .map2 <- suppressMessages(tryCatch(AnnotationDbi::mapIds(
        .odb, keys = .todo, keytype = "ALIAS", column = "ENTREZID",
        multiVals = "first"), error = function(e) setNames(rep(NA_character_, 
                                                               length(.todo)), 
                                                           .todo)))
      .map1[.todo] <- .map2[.todo]
    }
    sym2eg <- .map1[!is.na(.map1)]
    sym2eg <- sym2eg[!duplicated(sym2eg)]          # one symbol per Entrez
    UNIVERSE <- unname(sym2eg)                     # universe: genes after detection
    .rate <- length(sym2eg) / length(genes)
    cat("\n--- 8.1 identifier mapping ---\n")
    cat("tested genes:", length(genes), "| mapped to Entrez:", length(sym2eg),
        sprintf("(%.0f%%)", 100 * .rate), "\n")
    if (.rate < 0.5)
      cat("  WARNING: under half of your genes mapped. Check that",
          GO_ORGDB, "is the right organism and that rownames(counts)\n",
          "  are gene SYMBOLS and not Ensembl/RefSeq IDs. GO results\n",
          "  below are not trustworthy until this is fixed.\n")
    
    ## 
    ## 8.2 GO gene sets ----
    ## 
    ## GO2ALLEGS, not GO2EG: the ALL form propagates genes up the
    ## ontology to ancestor terms, which is what makes a parent term
    ## testable. Sets are then intersected with the universe and
    ## sized after that intersection -- sizing before it inflates
    ## every p-value.
    
    .g2a  <- get(paste0(sub("\\.db$", "", GO_ORGDB), "GO2ALLEGS"))
    .sets <- AnnotationDbi::as.list(.g2a)
    .sets <- lapply(.sets, function(x) unique(x[!is.na(x)]))
    
    .ann <- suppressMessages(AnnotationDbi::select(
      GO.db, keys = names(.sets), 
      keytype = "GOID", 
      columns = c("TERM", "ONTOLOGY")))
    .ann <- .ann[!duplicated(.ann$GOID), ]
    .keep_ont <- .ann$GOID[.ann$ONTOLOGY %in% GO_ONT]
    .sets <- .sets[intersect(names(.sets), .keep_ont)]
    
    GO_SETS <- lapply(.sets, function(x) intersect(x, UNIVERSE))
    .sz <- lengths(GO_SETS)
    GO_SETS <- GO_SETS[.sz >= GO_MINSIZE & .sz <= GO_MAXSIZE]
    GO_TERM <- setNames(.ann$TERM, .ann$GOID)
    
    cat("\n--- 8.2 gene sets ---\n")
    cat("ontology:", paste(GO_ONT, collapse = "/"),
        "| terms after universe intersection and size filter:", 
        length(GO_SETS), "\n")
    cat("set size: median", median(lengths(GO_SETS)),
        "| range", paste(range(lengths(GO_SETS)), collapse = "-"), "\n")
    if (!length(GO_SETS)) stop("no GO sets survived filtering; loosen GO_MINSIZE/GO_MAXSIZE")
    
    ## 
    ## 8.3 the ranking statistic ----
    ## 
    ## One vector per region, over the same contrasts the selection
    ## used, so GO results and the panel are answering questions
    ## about the same comparison.
    
    .rank_tab <- do.call(rbind, lapply(GROUPS, function(g) {
      rv <- setdiff(GROUPS, g)
      lf <- LFC[, g, rv, drop = FALSE][, 1, ]
      if (is.null(dim(lf))) lf <- matrix(lf, ncol = 1, 
                                         dimnames = list(genes, rv))
      data.frame(gene = genes, region = g,
                 mean_lfc = rowMeans(lf, na.rm = TRUE),
                 min_lfc  = apply(lf, 1, min, na.rm = TRUE),
                 row.names = NULL, stringsAsFactors = FALSE)
    }))
    
    cat("\n--- 8.3 ranking statistic ---\n")
    cat("why mean_lfc is the default: a worst-case min over rivals is\n")
    cat("shifted down, so its negative tail is not interpretable.\n\n")
    print(as.data.frame(.rank_tab %>% 
                          group_by(region) %>%
                          summarise(mean_lfc_median = round(median(mean_lfc), 3),
                                    mean_lfc_pct_neg = round(100 * mean(mean_lfc < 0)),
                                    min_lfc_median  = round(median(min_lfc), 3),
                                    min_lfc_pct_neg = round(100 * mean(min_lfc < 0)), 
                                    .groups = "drop")),
          row.names = FALSE)
    cat("\nusing GO_RANK =", GO_RANK, "\n")
    
## 
## 8.4 GSEA ----
## 
## NES = ES / mean(|ES| of the permutation null of the same
## sign). It is the enrichment score corrected for gene-set size,
## which is what makes NES comparable across terms and across
## regions -- ES is not. fgseaMultilevel estimates the small
## p-values by adaptive multilevel sampling rather than by
## brute-force permutation, so p is not floored at 1/GO_NPERM.
##
## Sign convention: NES > 0 means the term sits at the top of
## that region's ranking, i.e. UP in this region relative to the
## other selection regions.

## One helper, used by both GSEA and the running-ES figure, so the
## two cannot drift apart.
##

.rank_vec <- function(st, g) {
  v <- st[[GO_RANK]]; s <- st$gene
  ok <- s %in% names(sym2eg) & is.finite(v)
  v <- v[ok]; eg <- unname(sym2eg[s[ok]])
  if (anyDuplicated(eg)) {
    o  <- order(abs(v), decreasing = TRUE)     # keep largest |stat|
    v  <- v[o]; eg <- eg[o]
    k  <- !duplicated(eg); v <- v[k]; eg <- eg[k]
  }
  sort(setNames(v, eg), decreasing = TRUE)
}
    
 set.seed(1)
 gsea <- do.call(rbind, lapply(GROUPS, function(g) {
   st <- .rank_tab[.rank_tab$region == g, ]
   v  <- .rank_vec(st, g)
   r  <- fgsea::fgsea(pathways = GO_SETS, stats = v,
                      minSize = GO_MINSIZE, maxSize = GO_MAXSIZE,
                      nPermSimple = GO_NPERM)
   r  <- as.data.frame(r)
   if (!nrow(r)) return(NULL)
   ## leading edge back to symbols, for a readable table
   eg2sym <- setNames(names(sym2eg), sym2eg)
   r$leading_genes <- vapply(r$leadingEdge, function(e)
     paste(head(eg2sym[e], 12), collapse = "/"), character(1))
   r$n_leading <- lengths(r$leadingEdge)
   data.frame(region = g, GOID = r$pathway, term = unname(GO_TERM[r$pathway]),
              NES = r$NES, ES = r$ES, pval = r$pval, padj = r$padj,
              size = r$size, n_leading = r$n_leading,
              leading_genes = r$leading_genes,
              row.names = NULL, stringsAsFactors = FALSE)
 }))
 gsea$leadingEdge_list <- NULL
 gsea <- gsea[order(gsea$region, gsea$padj, -abs(gsea$NES)), ]
 
 cat("\n--- 8.4 GSEA ---\n")
 cat("terms tested per region:", length(GO_SETS),
     "| NES sign: + = up in that region vs the other selection regions\n\n")
 print(as.data.frame(gsea %>% group_by(region) %>%
                       summarise(sig = sum(padj < GO_PADJ, na.rm = TRUE),
                                 sig_up = sum(padj < GO_PADJ & NES > 0, 
                                              na.rm = TRUE),
                                 sig_dn = sum(padj < GO_PADJ & NES < 0, 
                                              na.rm = TRUE),
                                 top_NES = round(max(NES, na.rm = TRUE), 2), 
                                 .groups = "drop")),
       row.names = FALSE)
 
 cat("\ntop", GO_TOP_PLOT, "terms per region by adjusted p:\n")
 print(as.data.frame(gsea %>% filter(padj < GO_PADJ, NES > 0) %>%
                       group_by(region) %>% 
                       slice_head(n = GO_TOP_PLOT) %>% 
                       ungroup() %>%
                       mutate(term = substr(term, 1, 46), 
                              NES = round(NES, 2),
                              padj = signif(padj, 2)) %>%
                       dplyr::select(region, term, NES, padj, size, n_leading)), 
       row.names = FALSE)
 
 ## 
 ## 8.5 redundancy pruning ----
 ## 
 ## GO is a DAG, so a real signal shows up as a dozen nested terms
 ## with near-identical leading edges. For the FIGURE only, keep
 ## the best-scoring term of each redundant cluster. The full
 ## table is exported unpruned.
 
 prune <- function(d, jac = GO_JACCARD) {
      d <- d[order(d$padj, -abs(d$NES)), ]
      if (nrow(d) < 2) return(d)
      le <- lapply(d$leading_genes, function(s) strsplit(s, "/", fixed = TRUE)[[1]])
      keep <- logical(nrow(d))
      for (i in seq_len(nrow(d))) {
        red <- FALSE
        for (j in which(keep)) {
          u <- length(union(le[[i]], le[[j]]))
          if (u > 0 && length(intersect(le[[i]], le[[j]])) / u >= jac) { red <- TRUE; break }
        }
        keep[i] <- !red
      }
      d[keep, ]
    }
 gsea_pruned <- do.call(rbind, lapply(split(gsea, gsea$region), prune))
 cat("\n--- 8.5 redundancy ---\n")
 cat("terms before pruning:", nrow(gsea), "| after (figure only):",
     nrow(gsea_pruned), "\n")
 
 ## 
 ## 8.6 ORA on the selected panel----
 ## 
 ## Hypergeometric, computed directly so the universe is explicit.
 ## Foreground = the genes actually in the panel for that region.
 ## Universe  = every detected, mapped gene (not all of GO, and
 ## NOT the whole genome -- both inflate significance).
 
 ora <- do.call(rbind, lapply(GROUPS, function(g) {
      fg <- panel$gene[panel$region == g]
      fg <- unname(sym2eg[intersect(fg, names(sym2eg))])
      if (length(fg) < 3) return(NULL)
      N <- length(UNIVERSE); K <- length(fg)
      res <- do.call(rbind, lapply(names(GO_SETS), function(id) {
        s <- GO_SETS[[id]]; M <- length(s); k <- length(intersect(fg, s))
        if (k < 2) return(NULL)
        data.frame(GOID = id, k = k, K = K, M = M, N = N,
                   fold = (k / K) / (M / N),
                   pval = stats::phyper(k - 1, M, N - M, K, lower.tail = FALSE),
                   stringsAsFactors = FALSE)
      }))
      if (is.null(res) || !nrow(res)) return(NULL)
      res$padj <- stats::p.adjust(res$pval, "BH")
      res$region <- g; res$term <- unname(GO_TERM[res$GOID])
      res$genes <- vapply(res$GOID, function(id) {
        eg2sym <- setNames(names(sym2eg), sym2eg)
        paste(head(eg2sym[intersect(fg, GO_SETS[[id]])], 12), collapse = "/")
      }, character(1))
      res[order(res$padj), c("region", "GOID", "term", "k", "K", "M", "N",
                             "fold", "pval", "padj", "genes")]
    }))
 
 cat("\n--- 8.6 ORA on the panel genes ---\n")
 if (is.null(ora) || !nrow(ora)) {
   cat("no region had >=3 mapped panel genes; ORA skipped.\n")
 } else {
   cat("universe:", length(UNIVERSE), "detected+mapped genes",
       "| multiple testing: BH within region\n\n")
   print(as.data.frame(ora %>% group_by(region) %>%
                         summarise(panel_mapped = dplyr::first(K),
                                   sig = sum(padj < GO_PADJ), .groups = "drop")), row.names = FALSE)
   cat("\ntop ORA terms per region:\n")
   print(as.data.frame(ora %>% filter(padj < GO_PADJ) %>% group_by(region) %>%
                         slice_head(n = 3) %>% ungroup() %>%
                         mutate(term = substr(term, 1, 42), fold = round(fold, 1),
                                padj = signif(padj, 2)) %>%
                         dplyr::select(region, term, k, M, fold, padj)), row.names = FALSE)
   
   ## do the two methods agree?
   cmp <- merge(ora[, c("region", "GOID", "padj")],
                gsea[, c("region", "GOID", "padj", "NES")],
                by = c("region", "GOID"), suffixes = c("_ora", "_gsea"))
   if (nrow(cmp)) {
     cat("\nagreement (terms significant by each method):\n")
     print(as.data.frame(cmp %>% group_by(region) %>% summarise(
       both = sum(padj_ora < GO_PADJ & padj_gsea < GO_PADJ),
       ora_only = sum(padj_ora < GO_PADJ & padj_gsea >= GO_PADJ),
       gsea_only = sum(padj_ora >= GO_PADJ & padj_gsea < GO_PADJ),
       .groups = "drop")), row.names = FALSE)
     cat("  ora_only terms rest on a few top genes; treat with care.\n")
   }
 }
   
##
## 8.7 figures ----
##
## (a) NES heatmap: union of each region's top terms, every
##     region's NES shown for all of them, so a term that is
##     specific reads as one coloured cell in its row.
## (b) per-region lollipop of NES with the leading-edge count.
 
.sel_terms <- gsea_pruned %>% 
   filter(padj < GO_PADJ, NES > 0) %>% 
   arrange(-NES) %>% 
   group_by(region) %>% 
   slice_head(n = GO_TOP_PLOT) %>% 
   ungroup()
 
 .tlr_terms <- gsea_pruned |>
   filter(grepl("^toll-like receptor \\d+ signaling pathway", term)) |>
   group_by(region) |>
   summarise(
     n_TLR_pathways = n(),
     mean_NES = mean(NES, na.rm = TRUE),
     median_NES = median(NES, na.rm = TRUE),
     min_FDR = min(padj, na.rm = TRUE),
     n_significant = sum(padj < 0.05, na.rm = TRUE)
   )
 
 if (nrow(.sel_terms)) {
   hm <- gsea %>% filter(GOID %in% .sel_terms$GOID) %>%
     mutate(region = factor(region, levels = GROUPS),
            label  = paste0(substr(term, 1, 44),
                            ifelse(nchar(term) > 44, "\u2026", "")))
   ord <- .sel_terms %>% 
     arrange(match(region, GROUPS), desc(NES))
   hm$label <- factor(hm$label,
                      levels = rev(unique(paste0(substr(ord$term, 1, 44),
                                                 ifelse(nchar(ord$term) > 44, 
                                                        "\u2026", "")))))
   hm$star <- ifelse(is.na(hm$padj), "", ifelse(hm$padj < 0.001, "***",
                                                ifelse(hm$padj < 0.01, "**", 
                                                       ifelse(hm$padj < GO_PADJ, 
                                                              "*", 
                                                              ""))))
   .lim <- max(abs(hm$NES), na.rm = TRUE)
   
   
   
   .tlr_terms$star<- ifelse(is.na(.tlr_terms$min_FDR), "", ifelse(.tlr_terms$min_FDR < 0.001, "***",
                                                     ifelse(.tlr_terms$min_FDR < 0.01, "**", 
                                                            ifelse(.tlr_terms$min_FDR < GO_PADJ, 
                                                                   "*", 
                                                                   ""))))
   
   
   tlr_terms <- as.data.frame(.tlr_terms) |> 
     mutate(NES = median_NES, 
            label = 'TLR signalling') |> 
     dplyr::select(region, NES, star, label)
   
   ord <- tlr_terms %>% 
     arrange(match(region, GROUPS), desc(NES))
   
   rbind(hm |> dplyr::select(region, NES, star, label),
         tlr_terms |> dplyr::select(region, NES, star, label)) |>
     arrange(region, desc(NES)) -> df_plot
   
   p_nes <- df_plot |>
     mutate(label = factor(label, levels = unique(df_plot$label)))|> 
     ggplot(aes(region, label, fill = NES)) +
     geom_tile(colour = "white", linewidth = .4) +
     geom_text(aes(label = ifelse(star == "", 
                                  sprintf("%.1f", NES),
                                  paste0(sprintf("%.1f", NES), "\n", 
                                         star))),
               size = 2.5, lineheight = .8,
               colour = ifelse(abs(df_plot$NES) > .7 * .lim, "white", "grey15")) +
     scale_fill_gradient2(high = "darkblue", mid = "grey96", low = "#C8A21C",
                          midpoint = 0, limits = c(-.lim, .lim),
                          name = "NES") +
     labs(x = NULL, y = NULL) +
     theme_minimal(base_size = 9) +
     theme(panel.grid = element_blank(),
           axis.text.y = element_text(size = 7.5),
           axis.text.x = element_text(face = "bold"),
           plot.title.position = "plot",
           plot.title = element_text(size = 10),
           plot.subtitle = element_text(size = 7.5, colour = "grey30"))
   

   .h    <- max(3.2, 0.26 * nlevels(hm$label) + 1.9)
   .lab  <- max(nchar(levels(hm$label)))
   .wlab <- 0.055 * .lab                       
   
   ggsave(file.path(OUT_FIG, "go_nes_heatmap.png"), 
          plot = p_nes,
          width = max(5.2, .wlab + 0.55 * length(GROUPS) + 1.5),
          height = .h, dpi = 300)
   svglite::svglite(paste0(OUT_FIG,'/GO_NES.svg'),
                    width = max(5.2, .wlab + 0.55 * length(GROUPS) + 1.5),
                    height = .h)
   print(p_nes)
   dev.off()
   
   cat("\nwrote go_nes_heatmap.png (", nlevels(hm$label), "terms x",
       length(GROUPS), "regions )\n")
   
   lp <- .sel_terms %>% mutate(region = factor(region, levels = GROUPS),
                               label = paste0(substr(term, 1, 40), ifelse(nchar(term) > 40, "\u2026", ""))) %>%
     arrange(region, NES) %>% mutate(row = row_number())

   p_lp <- ggplot(lp, aes(NES, reorder(paste(row), row))) +
     geom_segment(aes(x = 0, xend = NES, yend = reorder(paste(row), row)),
                  colour = "grey75", linewidth = .4) +
     geom_point(aes(size = n_leading, colour = -log10(padj))) +
     facet_wrap(~region, scales = "free_y", nrow = 3) +
     scale_y_discrete(breaks = paste(lp$row), labels = lp$label) +
     scale_colour_viridis_c(option = "mako", direction = -1,
                            name = expression(-log[10] * " BH p")) +
     scale_size_continuous(range = c(1.2, 3.6), name = "leading-edge genes") +
     expand_limits(x = 0.1) +
     labs( x = "NES", y = NULL) +
     theme_minimal(base_size = 9) +
     theme(panel.grid.major.y = element_blank(),
           axis.text.y = element_text(size = 6.8),
           strip.text = element_text(face = "bold"),
           plot.subtitle = element_text(size = 7.5, colour = "grey30"))
   ggsave(file.path(OUT_FIG, "go_nes_dotplot.png"), plot = p_lp,
          width = 4.0 + 0.048 * max(nchar(lp$label)) * 2,
          height = 1.5 + 1.05 * ceiling(length(GROUPS) / 2) * GO_TOP_PLOT / 2.4,
          dpi = 300)
   cat("wrote go_nes_dotplot.png\n")
   

   best <- .sel_terms %>% group_by(region) %>% arrange(-NES) %>% slice_head(n = 1) %>% ungroup()
   cur <- do.call(rbind, lapply(seq_len(nrow(best)), function(i) {
     g  <- best$region[i]
     st <- .rank_tab[.rank_tab$region == g, ]
     v  <- .rank_vec(st, g)
     hit <- names(v) %in% GO_SETS[[best$GOID[i]]]
     nh <- sum(hit); nm <- length(v) - nh
     if (!nh) return(NULL)
     inc <- ifelse(hit, abs(v)^1 / sum(abs(v[hit])), -1 / nm)
     data.frame(region = g, rank = seq_along(v), run = cumsum(inc),
                hit = hit, term = substr(best$term[i], 1, 34),
                nes = best$NES[i], stringsAsFactors = FALSE)
   }))
   if (!is.null(cur) && nrow(cur)) {
     cur$region <- factor(cur$region, levels = GROUPS)
     lab <- cur %>% group_by(region) %>%
       summarise(term = dplyr::first(term), nes = dplyr::first(nes),
                 y = max(run), .groups = "drop")
     p_run <- ggplot(cur, aes(rank, run)) +
       geom_hline(yintercept = 0, colour = "grey80", linewidth = .3) +
       geom_line(colour = "#B2182B", linewidth = .5) +
       geom_rug(data = subset(cur, hit), sides = "b",
                length = unit(.045, "npc"), colour = "grey35", linewidth = .2) +
       geom_text(data = lab, aes(x = Inf, y = y, label =
                                   paste0(term, "\nNES ", sprintf("%.2f", nes))),
                 hjust = 1.02, vjust = 1, size = 2.4, colour = "grey20") +
       facet_wrap(~region, nrow = 3) +
       labs(x = "Gene rank", y = "Running ES") +
       theme_classic(base_size = 9) +
       theme(panel.grid.minor = element_blank(),
             strip.text = element_text(face = "bold"),
             plot.subtitle = element_text(size = 7.5, colour = "grey30"))
     ggsave(file.path(OUT_FIG, "go_running_es.png"), plot = p_run,
            width = 8.4, height = 2.1 * ceiling(length(GROUPS) / 3) + .9, dpi = 300)
     
     svglite::svglite(file.path(OUT_FIG, 'go_running_es.svg'),
                       height = 2.1 * ceiling(length(GROUPS) / 3) + .9,
                       width = 8.4)
     print(p_run)
     dev.off()
     
     cat("wrote go_running_es.png\n")
   }
 } else {
   cat("\nno term reached BH p <", GO_PADJ, "in any region; no GO figure written.\n")
 }
 
 ## 
 ## 8.8 exports ----
 ## 
 write.csv(gsea, file.path(OUT_CSV, "go_gsea_all.csv"), row.names = FALSE)
 nes_wide <- gsea %>% dplyr::select(region, GOID, term, NES) %>%
   tidyr::pivot_wider(names_from = region, values_from = NES)
 write.csv(nes_wide, file.path(OUT_CSV, "go_nes_matrix.csv"), row.names = FALSE)
 if (!is.null(ora) && nrow(ora))
   write.csv(ora, file.path(OUT_CSV, "go_ora_panel.csv"), row.names = FALSE)
 cat("\nwrote go_gsea_all.csv, go_nes_matrix.csv",
     ifelse(!is.null(ora) && nrow(ora), ", go_ora_panel.csv", ""), "\n", sep = "")
 
 cat("\n--- how to read section 8 ---\n")
 cat("* NES is comparable across terms and regions; ES is not.\n")
 cat("* NES sign is relative to the OTHER selection regions, not\n")
 cat("  to the tumour controls -- those were never in the contrast.\n")
 cat("* a term significant in two regions is not a contradiction:\n")
 cat("  GSEA asks about a shift, and two regions can both shift up\n")
 cat("  relative to the remaining ones.\n")
 cat("* ORA inherits LFC_MIN and N_TOP. If you change those, ORA\n")
 cat("  moves and GSEA barely does. That asymmetry is the point.\n")
 
  }   
}
#Fin
