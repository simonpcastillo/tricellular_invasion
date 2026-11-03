# Load the repository configuration; set PVN_CONFIG to use a different file.
.bootstrap_path <- grep('^--file=', commandArgs(FALSE), value = TRUE)
if (!length(.bootstrap_path)) stop('Run with Rscript, not source(): cannot infer script path')
.bootstrap_path <- normalizePath(sub('^--file=', '', .bootstrap_path[1]), mustWork = TRUE)
.project_root <- normalizePath(file.path(dirname(.bootstrap_path), '..'), mustWork = TRUE)
.config_file <- Sys.getenv('PVN_CONFIG', file.path(.project_root, 'config', 'config.R'))
source(.config_file, local = FALSE)

################################################################################
# 2_1_pvn_final.R
# Project: Perivascular invasion paper
# Author: Simon Castillo
# Last update: 17/09/2026
################################################################################
##
##
# 0. preamble ----
##

dir.create(OUT_LOG, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(OUT_LOG, '01_pvn_discovery.log')
# Fail before loading large datasets if the optional project helper code is absent.
.required_helpers <- c('pvn_enrichment_tests.R', 'pvn_vs_global_tests.R',
                       'pvn_composition_heatmap.R')
.missing_helpers <- .required_helpers[!file.exists(file.path(HELPERS_DIR, .required_helpers))]
if (length(.missing_helpers)) stop(
  'Missing required helper script(s): ', paste(.missing_helpers, collapse = ', '),
  '\nPlace the original files in ', HELPERS_DIR)
for (.f in FILES) if (!file.exists(.f)) stop('Missing QuPath input: ', .f)
sink(log_file, split = TRUE)
timestamp()
cat("Starting script: 2_1_pvn_final.R\n")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(RANN)
  library(spatstat.geom)
  library(scales)
  library(scattermore)
  library(svglite)
  library(ggh4x)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

library(cluster)
.ch_index <- function(X, lab) {
  k <- length(unique(lab)); n <- nrow(X)
  if (k < 2 || k >= n) return(NA_real_)
  tot <- sum(scale(X, scale = FALSE)^2)
  wss <- sum(vapply(split(seq_len(n), lab), function(ix){
    xs <- X[ix,,drop=FALSE]; sum(sweep(xs,2,colMeans(xs))^2)}, numeric(1)))
  ((tot - wss)/(k-1)) / (wss/(n-k))
}
separation_of <- function(idx, K) {                # per-image separation at K
  h  <- H[idx,,drop=FALSE]; pc <- prcomp(h, center=TRUE, scale.=FALSE)
  xp <- pc$x[, seq_len(min(MAX_PCS, ncol(pc$x))), drop=FALSE]
  d  <- dist(xp); lab <- cutree(hclust(d, method="ward.D2"), k=K)
  c(sil = mean(silhouette(lab, d)[,3]), ch = .ch_index(xp, lab))
}
  


# Input/output paths are set in config/config.R.
if (!dir.exists(OUT_CSV)) {
  dir.create(OUT_CSV, showWarnings = FALSE, recursive = TRUE)
}
if (!dir.exists(OUT_FIG)) {
  dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
}


set.seed(SEED)


## 
# 1. Load data ----
## 

df0 <- rbindlist(lapply(names(FILES), function(nm) {
  d <- fread(FILES[[nm]])
  d[, Image := nm][]
}), fill = TRUE)

setnames(df0, grep("^Centroid X", names(df0), value = TRUE), "cx")
setnames(df0, grep("^Centroid Y", names(df0), value = TRUE), "cy")

stopifnot(all(c("cx", "cy", "Classification") %in% names(df0)))
cat("\n== INPUT ==\n"); print(df0[, .N, by = Image])


##
# 2. Data tidyup ----
##

df0 <- df0[!Classification %in% c("Ependymal", "not_assigned", "multiple")]

df0[, class := fcase(
  Classification %in% c("IpEC", "mEC", "BEC"),  "Endothelial",
  Classification %in% NEFTEL, if (SPLIT_NEFTEL) Classification else "Tumour",
  Classification %in% c("OPC", "Oligodendro"),   "Oligodendro_OPC",
  default = Classification
)]
df0[, gid := .I]

cat("\n== CLASSES BY IMAGE ==\n")
print(dcast(df0[, .N, by = .(Image, class)], 
            class ~ Image, 
            value.var = "N", 
            fill = 0))

##
# 3. Region from signed distances ----
##

dist_cols <- grep("^Signed distance to annotation", 
                  names(df0), 
                  value = TRUE)

region_of <- function(d) {
  if (!length(dist_cols)) return(rep("unassigned", nrow(d)))
  M  <- as.matrix(d[, ..dist_cols])
  nm <- sub("^Signed distance to annotation with\\s*", "", dist_cols)
  nm <- tolower(trimws(gsub("[^A-Za-z ]+m?$", "", nm)))
  nm <- gsub("\\s+", "_", trimws(nm))
  
  inside <- M < 0 # negative signed distance = inside
  out <- rep("unassigned", nrow(M))
  for (r in c("core", "infiltrative_margin", "peritumoral")) {
    j <- which(nm == r)
    if (length(j)) out[out == "unassigned" & inside[, j[1]]] <- r
  }
  out
}
df0[, region := region_of(.SD), .SDcols = names(df0)]

if (mean(df0$region == "unassigned") > 0.5)
  warning("More than half of cells are region-unassigned")

cat("\n== REGION ==\n")
print(dcast(df0[, .N, by = .(Image, region)], 
            region ~ Image, 
            value.var = "N", 
            fill = 0))

# 
# 4. Neighbourhood composition ----
# 

NB_CLASSES <- sort(unique(df0[class != "Endothelial", class]))

build_composition <- function(d, radius, min_nb) {
  pl <- rbindlist(lapply(unique(d$Image), function(img) {
    s  <- d[Image == img]
    pp <- ppp(s$cx, s$cy, window = owin(range(s$cx), range(s$cy)), marks = s$class)
    cp <- closepairs(pp, r = radius, twice = FALSE)
    if (!length(cp$i)) return(NULL)
    I <- s$gid[cp$i]; J <- s$gid[cp$j]
    ec <- s[class == "Endothelial", gid]
    keep <- I %in% ec | J %in% ec
    I <- I[keep]; J <- J[keep]
    data.table(Image = img,
               e  = fifelse(I %in% ec, I, J),
               nb = fifelse(I %in% ec, J, I))
  }))
  if (!nrow(pl)) stop("no EC neighbour pairs at radius ", radius)
  pl[, nclass := d$class[nb]]
  pl <- pl[nclass != "Endothelial"] # EC define the anchor
  
  cn <- dcast(pl[, .N, by = .(Image, e, nclass)], 
              Image + e ~ nclass,
              value.var = "N", fill = 0)
  for (cc in setdiff(NB_CLASSES, names(cn))) cn[[cc]] <- 0L
  setcolorder(cn, c("Image", "e", NB_CLASSES))
  cn[, total := rowSums(.SD), .SDcols = NB_CLASSES]
  cn[total >= min_nb]
}

  cnt <- build_composition(df0, R_NB, MIN_NB)
  cnt <- merge(
    cnt,
    df0[class == "Endothelial", .(Image, e = gid, region)],
    by = c("Image", "e"),
    all.x = TRUE
  )
mat  <- as.matrix(cnt[, ..NB_CLASSES])
prop <- mat / rowSums(mat)

retention <- merge(
  df0[class == "Endothelial", .(n_ec_total = .N), by = Image],
  cnt[, .(n_ec_used = .N), by = Image], by = "Image")
retention[, pct_used := round(100 * n_ec_used / n_ec_total, 1)]

cat("\n== EC RETENTION ==\n")
print(retention)
cat(sprintf("Neighbourhood size: median %d, IQR %d-%d cells across %d classes\n",
            median(cnt$total), 
            quantile(cnt$total, .25), 
            quantile(cnt$total, .75),
            length(NB_CLASSES)))
cat(sprintf("Composition matrix is %.0f%% zeros; median %d of %d classes present\n",
            100 * mean(prop == 0), 
            median(rowSums(prop > 0)), 
            ncol(prop)))


# 
# 5. Hellinger transform ----
# 

hellinger <- function(p) sqrt(p)

H <- hellinger(prop)
H <- H[, apply(H, 2, sd) > 1e-8, drop = FALSE]

pca   <- prcomp(H, center = TRUE, scale. = FALSE)
X_pca <- pca$x[, seq_len(min(MAX_PCS, ncol(pca$x))), 
               drop = FALSE]

cat(sprintf("\nPCs retained: %d (%.1f%% of variance)\n", 
            ncol(X_pca),
            100 * sum(pca$sdev[seq_len(ncol(X_pca))]^2) / sum(pca$sdev^2)))


# 
# 6. Independent per image clustering ----
# 
# Here each image is clustered on its own and the centroids are matched later

cluster_one <- function(idx, K) {
  h  <- H[idx, , drop = FALSE]
  pc <- prcomp(h, center = TRUE, scale. = FALSE)
  xp <- pc$x[, seq_len(min(MAX_PCS, ncol(pc$x))), drop = FALSE]
  cutree(hclust(dist(xp), method = "ward.D2"), k = K)
}

centroids_of <- function(idx, lab) {
  t(apply(prop[idx, , drop = FALSE], 2, function(v) tapply(v, lab, mean)))
}

IMGS <- unique(cnt$Image)
stopifnot(length(IMGS) == 2)
i1 <- which(cnt$Image == IMGS[1]); i2 <- which(cnt$Image == IMGS[2])

## Shuffled-centroid null: how well do centroids match by chance?
## With few classes, any two centroid vectors correlate highly because both are
## dominated by the same abundant classes. This calibrates "matched".
centroid_null_q95 <- function(A, B, n = N_PERM_CENTROID) {
  quantile(replicate(n, max(apply(cor(A, B[sample(nrow(B)), ]), 1, max))),
           1 - REPRO_ALPHA)
}

K_eval <- rbindlist(lapply(K_CANDIDATES, function(K) {
  l1 <- cluster_one(i1, K); l2 <- cluster_one(i2, K)
  A <- centroids_of(i1, l1); B <- centroids_of(i2, l2)
  cm <- cor(A, B)
  best <- apply(cm, 1, max)
  q95 <- centroid_null_q95(A, B)
  sizes <- c(table(l1), table(l2))
  prop1 <- table(l1) / length(l1); prop2 <- table(l2) / length(l2)
  s1 <- separation_of(i1, K); s2 <- separation_of(i2, K)
  mean_sil = mean(c(s1["sil"], s2["sil"]))  
  min_sil = min(s1["sil"], s2["sil"])
  mean_CH  = mean(c(s1["ch"],  s2["ch"]))
  data.table(K = K,
             n_replicating = sum(best > q95),
             median_best_r = median(best),
             null_q95      = q95,
             min_size      = min(sizes),
             min_img_prop  = min(c(prop1, prop2)),
             mean_sil = mean_sil,
             min_sil=min_sil,
             mean_CH=mean_CH)
}))

cat("\n== K EVALUATION: replication against a shuffled-centroid null ==\n")
print(K_eval)

K_ok  <- K_eval[min_size >= MIN_CLUSTER_N & min_img_prop >= MIN_IMG_PROP_T & min_sil > 0.10]
if (!nrow(K_ok)) stop('No K candidate passes minimum size, proportion and silhouette thresholds.')
K_SEL <- K_ok[order(-mean_sil, -n_replicating),][1,K]
cat("\nSELECTED K =", K_SEL, "\n")


lab1 <- cluster_one(i1, K_SEL) 
lab2 <- cluster_one(i2, K_SEL)
C1 <- centroids_of(i1, lab1)   
C2 <- centroids_of(i2, lab2)
cormat <- cor(C1, C2)
null_q95_empirical <- centroid_null_q95(C1, C2)
null_q95 <- if (is.na(REPLICATION_MIN_R)) null_q95_empirical else REPLICATION_MIN_R
cat('Replication threshold: ', null_q95,
    ' (empirical shuffled-centroid q95 = ', null_q95_empirical, ')\n', sep = '')

## Greedy one-to-one matching.
match_tab <- data.table()
cm <- cormat
while (nrow(match_tab) < min(dim(cm)) && any(is.finite(cm))) {
  w <- which(cm == max(cm, na.rm = TRUE), arr.ind = TRUE)[1, ]
  match_tab <- rbind(match_tab, data.table(
    c1 = rownames(cm)[w[1]], c2 = colnames(cm)[w[2]], r = cm[w[1], w[2]]))
  cm[w[1], ] <- NA; cm[, w[2]] <- NA
}
match_tab[, replicating := r > null_q95]
match_tab[, PVN := paste0("PVN_", seq_len(.N))]

cat("\n== CENTROID MATCHING (null q95 =", round(null_q95, 3), ") ==\n")
print(match_tab)
cat(sprintf("%d of %d niches replicate.\n", sum(match_tab$replicating), nrow(match_tab)))
if (!any(match_tab$replicating))
  warning("NO niche replicates above the shuffled-centroid null. ",
          "Downstream results describe within-image structure only.")

## Assign every retained EC to its matched niche label
pvn <- rep(NA_character_, nrow(cnt))
pvn[i1] <- match_tab$PVN[match(as.character(lab1), match_tab$c1)]
pvn[i2] <- match_tab$PVN[match(as.character(lab2), match_tab$c2)]
cnt[, PVN := pvn]
cnt[, replicating := PVN %in% match_tab[replicating == TRUE, PVN]]


##
# 7. Niche composition ----
##

prof <- as.data.table(prop)[, `:=`(
  Image  = cnt$Image,
  PVN    = cnt$PVN,
  region = cnt$region
)]

profile_by_image <- prof[
  , lapply(.SD, mean),
  by = .(Image, PVN, region),
  .SDcols = NB_CLASSES
]

profile_final <- prof[
  , c(lapply(.SD, mean), .(n = .N)),
  by = .(PVN, region),
  .SDcols = NB_CLASSES
]

setorder(profile_final, PVN, region)


cat("\n== NICHE COMPOSITION ==\n")
print(profile_final[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])


##
# 8. Neftel-state fraction per niche ----
##

build_neftel_counts <- function(d, radius) {
  pl <- rbindlist(lapply(unique(d$Image), function(img) {
    s  <- d[Image == img]
    pp <- ppp(s$cx, s$cy, window = owin(range(s$cx), range(s$cy)), marks = s$class)
    cp <- closepairs(pp, r = radius, twice = FALSE)
    if (!length(cp$i)) return(NULL)
    I <- s$gid[cp$i]; J <- s$gid[cp$j]; ec <- s[class == "Endothelial", gid]
    keep <- I %in% ec | J %in% ec; I <- I[keep]; J <- J[keep]
    data.table(Image = img,
               e  = fifelse(I %in% ec, I, J),
               nb = fifelse(I %in% ec, J, I))
  }))
  pl[, ncls := d$Classification[nb]]          # RAW neighbour subtype (not collapsed)
  pl <- pl[ncls %in% NEFTEL]                   # tumour neighbours only
  cn <- dcast(pl[, .N, by = .(Image, e, ncls)], Image + e ~ ncls,
              value.var = "N", fill = 0)
  for (cc in setdiff(NEFTEL, names(cn))) cn[[cc]] <- 0L
  cn[]
}

nef <- build_neftel_counts(df0, R_NB)

# align to the RETAINED EC set (same rows as cnt); ECs absent here had zero
# tumour neighbours -> fill their counts with 0
nef <- nef[cnt[, .(Image, e, PVN)], on = .(Image, e)]
for (cc in NEFTEL) nef[is.na(get(cc)), (cc) := 0L]
nef[, tum_tot := rowSums(.SD), .SDcols = NEFTEL]

# per-EC within-tumour fraction; NA when the EC has zero tumour neighbours
for (cc in NEFTEL)
  nef[, (paste0("neftel_", cc)) := fifelse(tum_tot > 0, get(cc) / tum_tot, NA_real_)]

fcols <- paste0("neftel_", NEFTEL)
neftel_profile <- nef[, c(lapply(.SD, mean, na.rm = TRUE),
                          .(n_ec_with_tumour = sum(tum_tot > 0))),
                      by = PVN, .SDcols = fcols]

# sanity checks (safe to keep; comment out once trusted
tum_expected <- if ("Tumour" %in% NB_CLASSES) cnt$Tumour else
  rowSums(as.matrix(cnt[, ..intersect(NEFTEL, NB_CLASSES)]))
#stopifnot(isTRUE(all.equal(nef$tum_tot, tum_expected)))                 # matches collapsed Tumour count
#stopifnot(isTRUE(all.equal(unname(rowSums(neftel_profile[, ..fcols])),  # fractions sum to 1 / niche
#                           rep(1, nrow(neftel_profile)))))

# ---- merge into profile_final ----------------------------------------------
profile_final <- merge(profile_final, neftel_profile, by = "PVN", all.x = TRUE)
setorder(profile_final, PVN)

cat("\n== NICHE COMPOSITION + within-tumour Neftel fractions ==\n")
print(profile_final[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])


profile_final %>%
  mutate(PVN = factor(str_remove_all(PVN, "_"), 
                      levels = paste0('PVN', 1:length(unique(profile_final$PVN))))) |> 
  dplyr::select(PVN, neftel_MES_like:neftel_NPC_like) %>% 
  pivot_longer(cols = neftel_MES_like:neftel_NPC_like, 
               names_to = 'neftel_state', 
               values_to = 'proportion') %>% 
  ggplot(aes(PVN, proportion, fill = neftel_state)) +
  geom_bar(stat = "identity") +
  theme_classic()-> stack_neftel_ec

svglite(file.path(OUT_FIG, 'pvn_stack_neftel_pvn.svg'), height = 2.5, width = 5.0)
print(stack_neftel_ec)
dev.off()

## 
# 9. Microglia enrichment against null ----
##
# The statistic is max(group mean) / global mean over K clusters fitted to
# maximise separation. It is bounded away from 1 by construction. The null
# reruns the whole pipeline with non-EC phenotype labels permuted within each
# image, preserving microglia counts, EC positions and total cell density.

MG <- if ("Microglia_core" %in% NB_CLASSES) "Microglia_core" else NA_character_

if (is.na(MG)) {
  warning("No Microglia_core class present; skipping the microglia test.")
  mg_obs <- NULL; mg_null <- NULL
} else {
  gl <- mean(prop[, MG])
  mg_obs <- prof[, .(n = .N, 
                     microglia = mean(get(MG))), 
                 by = PVN][
    , enrichment := microglia / gl][order(-microglia)]
  obs_max <- max(mg_obs$enrichment)
  
  pipeline_max <- function(d) {
    cn <- tryCatch(build_composition(d, R_NB, MIN_NB), 
                   error = function(e) NULL)
    if (is.null(cn) || nrow(cn) < 3 * K_SEL) return(NA_real_)
    m  <- as.matrix(cn[, ..NB_CLASSES])
    pr <- m / rowSums(m)
    if (mean(pr[, MG]) == 0) return(NA_real_)
    h  <- hellinger(pr)
    h <- h[, apply(h, 2, sd) > 1e-8, drop = FALSE]
    pc <- prcomp(h, center = TRUE, scale. = FALSE)
    xp <- pc$x[, seq_len(min(MAX_PCS, ncol(pc$x))), drop = FALSE]
    lb <- cutree(hclust(dist(xp), method = "ward.D2"), k = K_SEL)
    max(tapply(pr[, MG], lb, mean)) / mean(pr[, MG])
  }
  
  set.seed(SEED + 1)
  null_max <- vapply(seq_len(N_PERM_PIPELINE), function(i) {
    dd <- copy(df0)
    for (img in unique(dd$Image)) {
      ridx <- dd[Image == img & class != "Endothelial", which = TRUE]
      dd[ridx, class := sample(class)]
    }
    pipeline_max(dd)
  }, numeric(1))
  
  mg_null <- data.table(
    observed      = obs_max,
    null_median   = median(null_max, na.rm = TRUE),
    null_q95      = quantile(null_max, .95, na.rm = TRUE),
    ratio_to_null = obs_max / median(null_max, na.rm = TRUE),
    p_value       = (1 + sum(null_max >= obs_max, na.rm = TRUE)) / (N_PERM_PIPELINE + 1))
  
  cat("\n== MICROGLIA ENRICHMENT ==\n"); print(mg_obs)
  cat("\n== MICROGLIA vs PIPELINE NULL ==\n"); print(mg_null)
  cat(sprintf("\nReport this as %.2fx above the null median (raw ratio %.1fx is\n",
              mg_null$ratio_to_null, obs_max))
  cat("inflated by taking a maximum over clusters fitted to separate them).\n")
}

source(file.path(HELPERS_DIR, 'pvn_enrichment_tests.R'))
ec <- build_ec_table(cnt, prop, nef, df0)  # per-EC: PVN, region, mg_frac, neftel_*
# Preserve the original statistical tests; explicitly retain results for export.
mg_between_pvn <- NULL
mg_between_pvn_dual <- NULL
mg_between_pvn_glm <- NULL
neftel_group <- NULL
if (!is.na(MG)) {
  mg_between_pvn <- microglia_niche_enrichment(ec, n_perm = 2000)
  mg_between_pvn_dual <- run_dual(microglia_niche_enrichment, ec,
                                 n_perm = 2000, key = "PVN")
  mg_between_pvn_glm <- microglia_niche_glm(cnt, mg_class = MG)
}
neftel_between_pvn <- neftel_niche_enrichment(ec, strata = c("Image"), n_perm = 2000)
enriched <- if (!is.null(mg_between_pvn))
  mg_between_pvn[p_BH < 0.05 & effect_vs_rest > 0, PVN] else character(0)
if (length(enriched)) {
  neftel_group <- neftel_group_contrast(ec, enriched_pvns = enriched,
                                       strata = c("Image"), n_perm = 2000)
} else {
  message("No microglia-enriched PVNs at BH q < 0.05; skipping Neftel group contrast.")
}

source(file.path(HELPERS_DIR, 'pvn_vs_global_tests.R'))
pl <- build_pairs(df0, cnt, R_NB)  # EC-anchored non-EC neighbour pairs, retained ECs

# Microglia: perivascular fraction vs whole-tissue microglia fraction.
df_mg <- if (!is.na(MG))
  micro_vs_global(df0, cnt, pl, mg_class = MG, strata = c("Image")) else NULL

# Neftel: perivascular within-tumour state fraction vs whole-tumour composition.
df_nef <- neftel_vs_global(df0, cnt, pl, states = NEFTEL, strata = c("Image"))

##
# 10. Region first then subtype ----
##

cnt[, region  := df0$region[cnt$e]]
cnt[, subtype := df0$Classification[cnt$e]]
cnt[, cell_id := df0$cell_id[cnt$e]]


r2_of <- function(g) {
  g <- as.factor(g)
  Hs <- H
  gm <- colMeans(Hs)
  sst <- sum(sweep(Hs, 2, gm)^2)
  ssw <- sum(vapply(levels(g), function(l) {
    x <- Hs[g == l, , drop = FALSE]
    if (!nrow(x)) return(0)
    sum(sweep(x, 2, colMeans(x))^2)
  }, numeric(1)))
  1 - ssw / sst
}

variance_tab <- data.table(
  predictor = c("tumour region", 
                "annotated EC subtype", 
                "derived niche (ceiling)"),
  r2        = c(r2_of(cnt$region), 
                r2_of(cnt$subtype), 
                r2_of(cnt$PVN)))

cat("\n== VARIANCE IN NEIGHBOURHOOD COMPOSITION EXPLAINED ==\n")
print(variance_tab[, .(predictor, pct = round(100 * r2, 1))])
cat("(the derived niche is fitted to this variance.)\n")

region_tests <- rbindlist(lapply(IMGS, function(img) {
  rbindlist(lapply(unique(cnt[Image == img, region]), function(rg) {
    s <- cnt[Image == img & region == rg]
    tb <- table(s$PVN, s$subtype)
    if (nrow(tb) < 2 || ncol(tb) < 2 || sum(tb) < 20)
      return(data.table(Image = img, region = rg, n = sum(tb), p = NA_real_))
    data.table(Image = img, region = rg, n = sum(tb),
               p = suppressWarnings(chisq.test(tb, simulate.p.value = TRUE,
                                               B = 2000)$p.value))
  }))
}))

cat("\n== NICHE <-> EC SUBTYPE, WITHIN REGION ==\n")
print(region_tests)


##
# 11. Moran's I ----
##
## Normalisation
## Two-sided test centred on E[I] = -1/(n-1), 

moran_I <- function(x, edges, w) {
  n <- length(x); z <- x - mean(x)
  denom <- sum(z^2)
  if (!is.finite(denom) || denom <= 0 || nrow(edges) == 0) return(NA_real_)
  num <- 2 * sum(w * z[edges[, 1]] * z[edges[, 2]]) # each undirected edge twice
  (n / (2 * sum(w))) * num / denom          # W = 2*sum(w)
}

moran_test <- function(x, edges, w, nperm = N_PERM_MORAN) {
  obs <- moran_I(x, edges, w)
  if (!is.finite(obs)) return(list(I = NA_real_, p = NA_real_, EI = NA_real_))
  EI  <- -1 / (length(x) - 1)
  nul <- replicate(nperm, moran_I(sample(x), edges, w))
  p   <- (1 + sum(abs(nul - EI) >= abs(obs - EI), na.rm = TRUE)) / (nperm + 1)
  list(I = obs, p = p, EI = EI)
}

moran_res <- rbindlist(lapply(IMGS, function(img) {
  idx <- which(cnt$Image == img)
  co  <- as.matrix(df0[cnt$e[idx], .(cx, cy)])
  if (nrow(co) <= K_SPATIAL + 1) return(NULL)
  nn  <- nn2(co, k = K_SPATIAL + 1)
  ed  <- cbind(rep(seq_len(nrow(co)), K_SPATIAL), as.vector(nn$nn.idx[, -1]))
  dd  <- as.vector(nn$nn.dists[, -1])
  ed  <- t(apply(ed, 1, sort)); keep <- !duplicated(ed)
  ed  <- ed[keep, , drop = FALSE]; dd <- dd[keep]
  sig <- median(dd)
  w   <- exp(-dd^2 / (2 * sig^2))
  
  far <- dd > 2 * R_NB
  rbindlist(lapply(sort(unique(na.omit(cnt$PVN[idx]))), function(p) {
    x <- as.numeric(cnt$PVN[idx] == p)
    a <- moran_test(x, ed, w)
    b <- moran_test(x, ed[far, , drop = FALSE], w[far])
    data.table(Image = img, PVN = p, n = sum(x),
               I_all = a$I, p_all = a$p,
               I_nonoverlap = b$I, p_nonoverlap = b$p,
               pct_edges_overlapping = round(100 * mean(!far), 1))
  }))
}))

cat("\n== MORAN'S I (primary column: I_nonoverlap) ==\n")
print(moran_res)


##
# 12. Figures ----
##

theme_v2 <- theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        plot.title = element_text(face = "bold", size = 9))

heat_dt <- melt(profile_final[, c("PVN", NB_CLASSES), with = FALSE],
                id.vars = "PVN", variable.name = "class", 
                value.name = "proportion")

source(file.path(HELPERS_DIR, 'pvn_composition_heatmap.R'))
p_heatmap <- pvn_composition_heatmap(prop, cnt)
#p_heatmap <- pvn_composition_heatmap(prop_sep, cnt)


svglite(file.path(OUT_FIG, 'pvn_heatmap_scale.svg'), width = 3.0, height = 2.2)
print(p_heatmap)
dev.off()

p_repro <- ggplot(match_tab, 
                  aes(r, reorder(PVN, r), 
                      colour = replicating)) +
  geom_vline(xintercept = null_q95, 
             linetype = "22", 
             colour = "grey45") +
  geom_point(size = 3) +
  annotate("text", 
           x = null_q95, 
           y = 0.7, 
           label = " shuffled-centroid null (q95)",
           hjust = 0, size = 2.6, colour = "grey45") +
  scale_colour_manual(values = c(`TRUE` = "#2E7D32", `FALSE` = "#B23A48"),
                      name = "replicates") +
  scale_x_continuous(limits = c(min(0.3, min(match_tab$r)), 1.02)) +
  labs(title = "Cross-image replication of each niche",
       x = "Best-match centroid correlation", y = NULL) +
  theme_v2

p_microglia <- NULL
if (!is.null(mg_obs)) {
  p_microglia <- ggplot(mg_obs, aes(enrichment, 
                                    reorder(PVN, enrichment))) +
    geom_vline(xintercept = 1, colour = "grey70") +
    geom_vline(xintercept = mg_null$null_median, linetype = "22", colour = "#B23A48") +
    geom_col(fill = "grey70", width = 0.6) +
    labs(title = sprintf("Microglia enrichment: %.2fx above the null median (p = %.3f)",
                         mg_null$ratio_to_null, mg_null$p_value),
         x = "Microglia proportion / global mean", y = NULL) +
    theme_v2
}

var_dt <- copy(variance_tab)[, predictor := factor(predictor, levels = rev(predictor))]
p_variance <- ggplot(var_dt, aes(100 * r2, predictor)) +
  geom_col(fill = "grey70", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * r2)), hjust = -0.2, size = 2.8) +
  scale_x_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Variance in neighbourhood composition explained",
       x = "%", y = NULL) +
  theme_v2

plot_dt <- data.table(
  cx = df0$cx[cnt$e], cy = df0$cy[cnt$e],
  Image = cnt$Image, PVN = cnt$PVN, 
  replicating = cnt$replicating,
  dx_margin = df0$`Signed distance to annotation with tumour µm`[cnt$e],
  dx_core = df0$`Signed distance to annotation with core µm`[cnt$e],
  dx_peri = df0$`Signed distance to annotation with peritumoral µm`[cnt$e])

plot_dt |> 
  mutate(    PVN = factor(str_remove_all(PVN, '_'), 
                          levels = paste0("PVN", 1:length(unique(plot_dt$PVN))))
  ) -> plot_dt
setDT(plot_dt)  # dplyr mutate may strip data.table class

bg_dt <- df0[, .(cx, cy, Image)] 


p_map <- ggplot() +
  geom_point(data = bg_dt, 
             aes(cy / 1000, cx / 1000), 
             colour = "black", size = 0.15) +
  geom_point(data = plot_dt[!is.na(PVN)], aes(cy / 1000, cx / 1000, colour = PVN),
             size = 0.5) +
  coord_fixed() + scale_y_reverse() +
  facet_grid(Image ~ PVN, switch = "y") +
  scale_colour_viridis_d(option = "H", guide = "none") +
  labs(title = "", x = NULL, y = NULL) +
  theme_v2 + theme(axis.text = element_blank(), axis.ticks = element_blank())

p_summary <- (p_repro | p_variance) / p_heatmap
if (!is.null(p_microglia)) p_summary <- p_summary / p_microglia


### plot spatial maps ----
plt_mouse1 <-
  ggplot() +
  geom_scattermore(
    data = bg_dt %>% 
      filter(Image == "mouse1"), 
    aes(x = cy / 1000, y = cx / 1000),
    colour = "black",
    pointsize = 3,
    pixels = c(300, 300),
    interpolate = FALSE) +
  geom_point(
    data = plot_dt %>%
      filter(Image == "mouse1", !is.na(PVN)), 
    aes(cy / 1000, cx / 1000, colour = PVN), 
    size = 0.15) +
  coord_fixed() +
  scale_y_reverse() +
  facet_grid(. ~ PVN, switch = "y") +
  scale_colour_viridis_d(
    option = "H",
    guide = "none",
    begin = 0.1,
    end = 0.9) +
  labs( title = "", x = NULL, y = NULL) +
  theme_v2 +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
print(plt_mouse1)


svglite(file.path(OUT_FIG,"pvn_plt_distr_mouse1.svg"), width=5.5, height=1.2)
print(plt_mouse1)
dev.off()


plt_mouse2 <- ggplot() +
  geom_scattermore(
    data = bg_dt %>% filter(Image == "mouse2"), 
    aes( x = cy / 1000, y = cx / 1000),
    colour = "black",
    pointsize = 3,
    pixels = c(300, 300),
    interpolate = FALSE
    ) +
  geom_point(
    data = plot_dt %>%
      filter(Image == "mouse2", 
             !is.na(PVN)), 
    aes(cy / 1000, cx / 1000, 
        colour = PVN), 
    size = 0.15) +
  coord_fixed() +
  scale_x_reverse() +
  facet_grid(. ~ PVN, switch = "y") +
  scale_colour_viridis_d(
    option = "H",
    guide = "none",
    begin = 0.1,
    end = 0.9) +
  labs(title = "", x = NULL, y = NULL) +
  theme_v2 +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
print(plt_mouse2)


svglite(file.path(OUT_FIG, "pvn_plt_distr_mouse2.svg"), width = 5.5,height = 1.2)
print(plt_mouse2)
dev.off()


## plot density 2D ----
plot_dt %>% 
  mutate(
    PVN = factor(
      str_remove_all(PVN, '_'), 
      levels = paste0("PVN", 1:length(unique(plot_dt$PVN)))
    ),
    x = dx_peri / 1000,
    y = dx_core / 1000
  ) %>%
  ggplot(aes(x, y)) +
  geom_density2d(
    aes(
      colour = after_stat(
        case_when(
          x >= 0 & y >= 0 ~ "Infiltrative",
          x < 0  & y >= 0 ~ "PT",
          x >= 0 & y < 0  ~ "Bulk",
          TRUE ~ "NA"
        )
      )
    ), bins = 20,
    linewidth = 0.25
  ) +
  scale_colour_manual(
    values = c(
      "Infiltrative" = "#d2a6e2",
      "PT" = "gray60",
      "Bulk" = "#961212",
      "NA" = "gray80"
    )
  ) +
  facet_grid2(
    . ~ PVN,
    axes = "all"
    #scales = "free"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    linetype = 2
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    linetype = 2
  ) +
  theme_v2 +
  guides(color = "none") +
  labs(x = '', y = '') +
  theme(axis.text = element_blank())-> plt_density
print(plt_density)

svglite(file.path(OUT_FIG,'pvn_plt_density.svg'), width = 5.5, height = 1.2)
print(plt_density)
dev.off()

left_join(df_mg |> 
            mutate(effect_mg = effect,
                   p_BH_mg = p_BH) |> 
            dplyr::select(PVN, effect_mg, p_BH_mg),
          df_nef |> 
            #filter(state ==mes_state) |> 
            dplyr::select(PVN, state, effect, p_BH), by = 'PVN') |>
  mutate(
    PVN = factor(str_remove_all(PVN, '_'), 
                 levels = paste0("PVN", 1:length(unique(plot_dt$PVN)))),
    state_status = case_when(
      effect < 0 & p_BH < 0.05 ~ 'state_depleted',
      TRUE ~ 'stat_nondepleted',
    ),
    state = str_replace_all(state, "_", "-"))  |> 
    ggplot( aes(
      x = reorder(PVN, effect_mg),
      y = effect_mg,
      fill = state_status
    )) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
    scale_fill_manual(values = c('gray70', 'purple'))+
coord_flip() +
  theme_v2 +
  labs(
    x = NULL,
    y = "PVN vs WSI Microglia") +
    facet_wrap(.~state, axes = 'margins') +
    guides(fill = 'none') +
    theme(text = element_text(size = 6))-> plt_state

print(plt_state)

svglite(file.path(OUT_FIG,'pvn_plt_state_PVN.svg'), width = 2.4, height = 2.4)
print(plt_state)
dev.off()

##
# 13. EXPORT ----
##
  
fw <- function(x, nm) if (!is.null(x)) fwrite(x, file.path(OUT_CSV, nm))

fw(K_eval, "PVN_v2_K_evaluation.csv")
fw(match_tab, "PVN_v2_replication.csv")
fw(as.data.table(cormat, keep.rownames = "image1_cluster"),
 "PVN_v2_centroid_correlation.csv")
fw(profile_final,"PVN_v2_composition.csv")
fw(profile_by_image,"PVN_v2_composition_by_image.csv")
fw(mg_obs,"PVN_v2_microglia_enrichment.csv")
fw(mg_null,"PVN_v2_microglia_null.csv")
# Helper-generated enrichment tests (previously computed but not exported).
fw(mg_between_pvn, "PVN_v2_microglia_between_pvns.csv")
fw(mg_between_pvn_dual, "PVN_v2_microglia_marginal_region_adjusted.csv")
fw(mg_between_pvn_glm, "PVN_v2_microglia_glm.csv")
fw(neftel_between_pvn, "PVN_v2_neftel_between_pvns.csv")
fw(neftel_group, "PVN_v2_neftel_microglia_enriched_vs_rest.csv")
fw(df_mg, "PVN_v2_microglia_vs_whole_tissue.csv")
fw(df_nef, "PVN_v2_neftel_vs_whole_tissue.csv")
fw(variance_tab,"PVN_v2_variance_explained.csv")
fw(region_tests,"PVN_v2_niche_vs_subtype_within_region.csv")
fw(moran_res,"PVN_v2_morans_I.csv")
fw(retention,"PVN_v2_ec_retention.csv")
fw(cnt[, .(Image, cell = e, 
           cell_id = cell_id,
           PVN, 
           replicating, 
           region, 
           subtype, 
           total_neighbours = total)],
   "PVN_v2_assignments.csv")

ggsave(file.path(OUT_FIG, "PVN_v2_summary.pdf"),   p_summary, width = 11, height = 13)
ggsave(file.path(OUT_FIG, "PVN_v2_composition.pdf"), p_heatmap, width = 8, height = 5)
ggsave(file.path(OUT_FIG, "PVN_v2_maps.pdf"),      p_map,
       width = 3 + 2.2 * K_SEL, height = 6, limitsize = FALSE)

saveRDS(list(
  parameters = list(R_NB = R_NB, MIN_NB = MIN_NB, K_SEL = K_SEL, SEED = SEED,
                    transform = "hellinger", SPLIT_NEFTEL = SPLIT_NEFTEL),
  counts = cnt, prop = prop, hellinger = H, pca = pca,
  labels_image1 = lab1, labels_image2 = lab2,
  centroid_correlation = cormat, centroid_null_q95 = null_q95,
  match_table = match_tab, composition = profile_final,
  microglia = mg_obs, microglia_null = mg_null,
  variance = variance_tab, region_tests = region_tests,
  moran = moran_res, retention = retention
), file.path(OUT_CSV, "PVN_v2_model.rds"))

sink()
