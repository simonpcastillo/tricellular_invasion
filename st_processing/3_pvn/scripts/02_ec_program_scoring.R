# Load the repository configuration; set PVN_CONFIG to use a different file.
.bootstrap_path <- grep('^--file=', commandArgs(FALSE), value = TRUE)
if (!length(.bootstrap_path)) stop('Run with Rscript, not source(): cannot infer script path')
.bootstrap_path <- normalizePath(sub('^--file=', '', .bootstrap_path[1]), mustWork = TRUE)
.project_root <- normalizePath(file.path(dirname(.bootstrap_path), '..'), mustWork = TRUE)
.config_file <- Sys.getenv('PVN_CONFIG', file.path(.project_root, 'config', 'config.R'))
source(.config_file, local = FALSE)

################################################################################
# 02_ec_program_scoring.R
# Project: Perivascular invasion paper
# Author: Simon Castillo
# Last update: 22/09/2026
################################################################################

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix)
  library(ggplot2); library(patchwork)
  library(dplyr); library(stringr)
})


# Runtime settings, expression paths, reference GMT and output locations
# are defined in config/config.R (or PVN_CONFIG). Nothing is hardcoded here.

args <- commandArgs(trailingOnly = TRUE)
for (a in args) { kv <- strsplit(a, "=", fixed = TRUE)[[1]]
if (length(kv) == 2) assign(kv[1], type.convert(kv[2], as.is = TRUE), 
                            envir = .GlobalEnv) }
SEED <- EC_SEED
set.seed(SEED)
say <- function(...) cat(..., "\n", sep = "")
dir.create(OUT_RES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

find_col <- function(nms, patterns, what, required = TRUE) {
  for (p in patterns) { hit <- grep(p, nms, ignore.case = TRUE, value = TRUE)
  if (length(hit)) return(hit[1]) }
  if (required) stop(sprintf("Could not find the %s column. Available: %s",
                             what, paste(nms, collapse = ", ")), call. = FALSE)
  NA_character_
}

# 1. PVN assignments ----
if (!file.exists(PVN_ASSIGN)) stop('Missing PVN assignments: ', PVN_ASSIGN,
                                    '\nRun scripts/01_pvn_discovery.R first.')
assign_dt <- fread(PVN_ASSIGN)
stopifnot(all(c("Image", "cell_id", "PVN") %in% names(assign_dt)))
if (REPLICATING_ONLY && "replicating" %in% names(assign_dt)) {
  n0 <- nrow(assign_dt); assign_dt <- assign_dt[replicating == TRUE]
  say(sprintf("Anchors: %d of %d in replicating PVNs", nrow(assign_dt), n0))
}
assign_dt <- assign_dt[!PVN %in% c(NOT_ASSIGNED, EXCLUDE_PVN) & !is.na(PVN)]
assign_dt[, cid := as.character(cell_id)]
say("== EC anchors per PVN x Image ==")
print(dcast(assign_dt[, .N, by = .(Image, PVN)], 
            PVN ~ Image, value.var = "N", fill = 0))

# 2. load sparse expression, restrict to EC ----
load_image <- function(img) {
  if (!file.exists(EXPR[[img]])) stop('Missing sparse expression RDS: ', EXPR[[img]])
  obj  <- readRDS(EXPR[[img]])
  mat  <- obj$matrix
  if (is.null(mat) || is.null(obj$metadata))
    stop('Expected list(matrix=..., metadata=...) in ', EXPR[[img]])
  meta <- as.data.table(obj$metadata)
  id_m <- find_col(names(meta), 
                   c("^cell_id$", "^cell$", "^barcode$"), "cell id (expr)")
  meta[, cid := as.character(get(id_m))]
  meta <- meta[match(colnames(mat), cid)]
  stopifnot(identical(meta$cid, colnames(mat)))
  
  anc <- assign_dt[Image == img]
  keep_cid <- intersect(anc$cid, meta$cid)
  if (EC_CLASS_CHECK && LABEL_COL %in% names(meta)) {
    ec_cid   <- meta[get(LABEL_COL) == EC_CLASS, cid]
    n_pre    <- length(keep_cid)
    keep_cid <- intersect(keep_cid, ec_cid)
    if (length(keep_cid) < n_pre)
      message(sprintf("[%s] %d of %d anchors are not labelled %s in %s - dropped",
                      img, n_pre - length(keep_cid), n_pre, EC_CLASS, LABEL_COL))
  }
  lost <- sum(!anc$cid %in% meta$cid)
  if (lost) message(sprintf("[%s] %d PVN anchors have no expression profile - dropped",
                            img, lost))
  stopifnot(length(keep_cid) > 0)
  sub <- mat[, keep_cid, drop = FALSE]
  colnames(sub) <- paste(img, keep_cid, sep = "|")   # unique across images
  md  <- data.table(
    cell    = colnames(sub),
    .pvn    = anc$PVN[match(keep_cid, anc$cid)],
    .sec    = img,
    .reg    = if (REGION_COL %in% names(anc)) anc[[REGION_COL]][match(keep_cid, anc$cid)] else "all")
  md[is.na(.reg), .reg := "all"]
  say(sprintf("[%s] %d EC cells kept, %d genes", img, ncol(sub), nrow(sub)))
  list(mat = sub, md = md)
}

parts <- lapply(names(EXPR), load_image)
common_genes <- Reduce(intersect, lapply(parts, function(p) rownames(p$mat)))
stopifnot(length(common_genes) > 100)
M  <- do.call(cbind, lapply(parts, function(p) p$mat[common_genes, , drop = FALSE]))
md <- rbindlist(lapply(parts, `[[`, "md"))
setDF(md); rownames(md) <- md$cell
say(sprintf("Combined: %d EC cells x %d genes across %d sections",
            ncol(M), nrow(M), length(EXPR)))

# 3. build Seurat object with log-normalised input ----
so <- CreateSeuratObject(counts = M, 
                         assay = ASSAY, 
                         meta.data = md)

if (isTRUE(EXPR_IS_LOGNORM)) {
  SeuratObject::LayerData(so, 
                          layer = "data", 
                          assay = ASSAY) <- M
  say("EXPR_IS_LOGNORM: matrix placed in `data` layer, NormalizeData skipped")
} else {
  so <- NormalizeData(so, verbose = FALSE)
}
so$.pvn <- md$.pvn; so$.sec <- md$.sec; so$.reg <- md$.reg # assign pvn, smaple_id, and region
so$.stratum <- so$.sec   # section is the replication block
has_region <- uniqueN(so$.reg) > 1
if (length(unique(so$.pvn)) < 2) stop("need >= 2 niches after filtering")
genes <- rownames(so); lc <- tolower(genes)
NDET_COL <- paste0("nFeature_", ASSAY)

# OPT: curated EC programs ----

PROGRAMS <- list(
  arterial = c("Gkn3","Bmx","Efnb2","Sox17","Hey1","Vegfc","Cxcl12","Alpl","Sema3g",
               "Unc5b","Fbln5","Mgp","Clu","Htra1","Gja4","Gja5","Edn1","Cytl1",
               "Dkk2","Stmn2","Crip1","Tm4sf1","Ltbp4","Trim47","Nrp1"),
  capillary = c("Mfsd2a","Slc16a1","Rgcc","Car4","Tfrc","Slc7a5","Spock2","Cd300lg",
                "Slc6a6","Slco1a4","Abcb1a","Abcg2","Slc2a1","Aqp1","Kcne3","Itm2a"),
  venous = c("Nr2f2","Slc38a5","Vwf","Vcam1","Emcn","Lcn2","Ackr1","Il1r1","Plvap",
             "Tmem252","Lrg1","Icam1","Bgn","Adgrg6","Fmo2","Cfh","Cpe","Tsc22d1","Rbp7"),
  bbb_transport = c("Cldn5","Ocln","Tjp1","Tjp2","Marveld2","Cgn","Jam2","Jam3","Esam",
                    "Slc2a1","Slc7a5","Slc38a5","Slc16a1","Slc22a8","Slc6a6","Slco1a4",
                    "Slco2b1","Abcb1a","Abcg2","Abcc4","Mfsd2a","Lef1","Zic3","Foxq1",
                    "Tfrc","Car4","Slc39a10","Slc7a1"),
  tip_cell = c("Esm1","Apln","Angpt2","Cxcr4","Dll4","Kdr","Flt4","Nrp1","Pgf","Adm",
               "Lamb1","Marcksl1","Igfbp3","Ramp3","Nid2","Kcne3","Lxn","Trp53i11",
               "Plaur","Mcam","Robo4","Unc5b","Cd34","Jup","Ednrb"),
  stalk_quiescent = c("Aqp1","Cd36","Mgp","Sparcl1","Spock2","Jup","Klf2","Klf4","Id1",
                      "Sox18","Tek","Nos3","Thbd","Cav1","Ramp2","Cldn5","Bmp4","Fabp4"),
  proliferating = c("Mki67","Top2a","Ccnb1","Ccnb2","Birc5","Ube2c","Stmn1","Tubb5",
                    "Rrm2","Pclaf","Cdk1","Cdc20","Aurka","Plk1","Nusap1","Tpx2",
                    "Cenpf","Cenpe","Hmgb2","Tyms","Mcm3","Mcm5","Pcna","Ccna2",
                    "Cks1b","Racgap1","Hist1h2ap"),
  hypoxia_glycolysis = c("Vegfa","Slc2a1","Ldha","Pgk1","Aldoa","Hk2","Bnip3","Ddit4",
                         "Egln3","Ndrg1","Pdk1","Bnip3l","Adm","Ankrd37","P4ha1","P4ha2",
                         "Pfkl","Pfkp","Eno1","Tpi1","Slc16a3","Higd1a","Ero1a","Hilpda",
                         "Mif","Anxa2","Pgam1","Gpi1"),
  notch = c("Hey1","Hey2","Hes1","Heyl","Notch1","Notch4","Jag1","Jag2","Dll4",
            "Efnb2","Nrarp","Rbpj","Dtx1","Maml2"),
  ifn_antigen = c("B2m","H2-K1","H2-D1","H2-Q7","Cd74","H2-Ab1","H2-Aa","H2-Eb1",
                  "Tap1","Tap2","Psmb8","Psmb9","Isg15","Ifit1","Ifit2","Ifit3",
                  "Stat1","Stat2","Irf7","Irf9","Bst2","Oasl2","Rsad2","Usp18",
                  "Igtp","Iigp1","Zbp1","Gbp2","Gbp3","Cxcl10","Ifi27","Xaf1"),
  inflam_adhesion = c("Vcam1","Icam1","Sele","Selp","Ccl2","Cxcl1","Cxcl2","Nfkbia",
                      "Tnfaip3","Icam2","Ccl5","Cxcl9","Il6","Il1b","Traf1","Birc3",
                      "Nfkb1","Nfkb2","Relb","Sod2","Socs3","Ptgs2","Lcn2","Saa3","Cd40"),
  endmt = c("Acta2","Tagln","Cnn1","Fn1","Col1a1","Col1a2","Col3a1","Serpine1",
            "Snai1","Snai2","Twist1","Zeb1","Zeb2","Tgfb1","Tgfbr1","Vim","Thy1",
            "Postn","Sparc","Lox","Mmp2","Cdh2","S100a4","Tpm1"),
  ecm_basement = c("Col4a1","Col4a2","Lama4","Lamb1","Lamc1","Lama5","Nid1","Nid2",
                   "Hspg2","Fbln5","Emcn","Col18a1","Col15a1","Sparc","Timp3","Vwa1","Agrn"),
  hev = c("Chst4","Glycam1","Ccl21a","Ccl19","Fut7","Madcam1","Chst2","Gcnt1","St3gal4"),
  lymphatic = c("Prox1","Pdpn","Lyve1","Flt4","Mmrn1","Ccl21a","Nts","Reln","Stab1","Cldn11"),
  endothelial_adhesion = c("Vcam1","Icam1","Icam2","Sele","Selp","Cd40"),
  leukocyte_recruitment = c("Ccl2","Cxcl1","Cxcl2","Cxcl9","Cxcl10","Ccl5"),
  ifn_response = c("Isg15","Ifit1","Ifit2","Ifit3","Stat1","Stat2","Irf7","Irf9",
                   "Bst2","Oasl2","Rsad2","Usp18","Igtp","Iigp1","Zbp1","Gbp2",
                   "Gbp3","Cxcl10","Ifi27","Xaf1"),
  inflammatory_activation = c("Nfkb1","Nfkb2","Relb","Nfkbia","Tnfaip3","Traf1",
                              "Birc3","Il1r1","Il6","Il1b","Ptgs2","Socs3","Sod2",
                              "Nfkbiz","Nfkbia","Dusp1","Dusp4","Cxcl1","Cxcl2",
                              "Cxcl3","Ccl2","Lcn2","Saa3"),
  tlr_sensing = c("Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9",
                  "Myd88","Tirap","Ticam1","Ticam2","Irak1","Irak4","Traf6","TbK1"),
  tlr_activation = c("Nfkb1","Nfkb2","Rela","Nfkbia","Tnfaip3","Birc3","Socs3",
                     "Ptgs2","Il6","Cxcl1","Cxcl2","Cxcl3","Ccl2"),
  BEC = c('Ces2e','Piezo2','Ct55','Tc2n','Cenph','Cbln2','Disp2','Gm41452',
           'Tnfrsf9','Lincred1','Clca2','Troap','Cfap46','Cdkl1','Sdcbp2'),
  IpEC = c('Ccl24','Cd209f','Cd209a','Lyve1','Cd200r1','Tlr5','Mgl2','Clec10a',
            'Cd163','Rgs14','Gvin3','Gm36025','F13a1','Msr1','Cfp','C1qb'),
  mEC = c('Pln','Lancl3','Ptpru','Ceacam2','Ccdc69','Olfm5','Me3','Hbb-y',
           'Lurap1','Cytl1','Ssu2','Colca2','Ccn3','Slc27a3','C030037D09Rik',
           'Myh11','Zfp759','Sspn','Gm46851','Tcf15','Gm2102','Syne1','Syt7',
           'Crispld1','Gm15998','Gm38708','Odf3b','Npr2','Gm32918','Bmx')
)
if (!is.na(GMT) && file.exists(GMT)) {
  gl <- strsplit(readLines(GMT), "\t")
  nm_new <- vapply(gl, `[`, "", 1); n_rep <- sum(nm_new %in% names(PROGRAMS))
  for (g in gl) PROGRAMS[[g[1]]] <- setdiff(g[-(1:2)], "")
  say("Loaded ", length(gl), " program(s) from ", GMT, " (", n_rep,
      " replaced a built-in, ", length(gl) - n_rep, " new)")
}


# ==== match programs, build expression-matched random null banks ===========
cnt  <- SeuratObject::LayerData(so, layer = "counts", assay = ASSAY)
detf <- Matrix::rowSums(cnt > 0) / ncol(cnt)
map1 <- function(g) { i <- match(tolower(g), lc); genes[i[!is.na(i)]] }
P <- lapply(PROGRAMS, map1)
nm <- sapply(P, length) 
dropped <- names(P)[nm < MIN_GENES]
P <- P[nm >= MIN_GENES]
if (length(dropped)) say("Dropped (< ", MIN_GENES, " genes present): ",
                         paste(dropped, collapse = ", "))
if (!length(P)) stop("no program had >= ", MIN_GENES, " genes in this object")

set.seed(SEED)
so$.half <- "test"
hs <- as.data.table(so@meta.data)[, .(cell = colnames(so), .pvn, .stratum)]
dh <- hs[, .(cell = sample(cell, floor(.N / 2))), by = .(.pvn, .stratum)]$cell
so$.half[colnames(so) %in% dh] <- "derive"

if (isTRUE(DE_NOVO)) {
  sod <- so[, so$.half == "derive"]
  vf <- VariableFeatures(FindVariableFeatures(sod, 
                                              nfeatures = 2000, 
                                              verbose = FALSE))
  if (requireNamespace("RcppML", quietly = TRUE)) {
    X <- SeuratObject::LayerData(sod, layer = "data", assay = ASSAY)[vf, ]
    f <- RcppML::nmf(X, k = N_DE_NOVO, verbose = FALSE, seed = SEED)
    W <- f$w; rownames(W) <- vf
    for (k in seq_len(ncol(W)))
      P[[paste0("denovo_", k)]] <- rownames(W)[order(W[, k], decreasing = TRUE)][1:40]
    say("Added ", N_DE_NOVO, " de novo NMF programs (derived on held-in half)")
  } else {
    sod <- ScaleData(sod, features = vf, verbose = FALSE)
    sod <- RunPCA(sod, features = vf, npcs = N_DE_NOVO, verbose = FALSE)
    ld <- Loadings(sod, "pca")
    for (k in seq_len(ncol(ld))) {
      P[[paste0("denovo_PC", k, "_hi")]] <- rownames(ld)[order(ld[, k], decreasing = TRUE)][1:40]
      P[[paste0("denovo_PC", k, "_lo")]] <- rownames(ld)[order(ld[, k])][1:40]
    }
    say("RcppML absent; added ", ncol(ld), " de novo PCA programs")
  }
}
say("Programs scored: ", length(P))

szbin <- function(n) pmax(5L, as.integer(round(n / 5) * 5))
bins <- sort(unique(szbin(sapply(P, length)))); bank <- list()
for (bn in bins) {
  ref <- unlist(P[szbin(sapply(P, length)) == bn])
  lo <- quantile(detf[ref], 0.10, na.rm = TRUE); hi <- quantile(detf[ref], 0.90, na.rm = TRUE)
  poolg <- setdiff(names(detf)[detf >= lo & detf <= hi], unlist(P))
  if (length(poolg) < bn * 3) poolg <- setdiff(names(detf)[detf > 0], unlist(P))
  bank[[as.character(bn)]] <- replicate(NR, sample(poolg, bn), simplify = FALSE)
}
allsets <- c(P, unlist(lapply(names(bank), function(s)
  setNames(bank[[s]], paste0("__null_", s, "_", seq_along(bank[[s]])))), recursive = FALSE))

# ==== score ================================================================
say("Scoring ", length(allsets), " gene sets (", length(P), " real + ",
    length(allsets) - length(P), " matched null) ...")
so <- AddModuleScore(so, features = allsets, name = "MS", ctrl = CTRL, verbose = FALSE)
S <- as.data.table(so@meta.data)[, c(".pvn", ".stratum", ".reg", ".half",
                                     paste0("MS", seq_along(allsets))), with = FALSE]
setnames(S, paste0("MS", seq_along(allsets)), names(allsets))

L <- melt(S, id.vars = c(".pvn", ".stratum", ".reg", ".half"),
          variable.name = "set", value.name = "s")
is_dn <- function(x) grepl("^denovo_", as.character(x))

niche_eff <- function(Lx) {
  gr <- Lx[, .(m = mean(s), n = .N), by = .(set, .stratum, .pvn)][n >= MIN_CELLS]
  gr[, eff := m - mean(m), by = .(set, .stratum)]
  E <- gr[, .(eff = mean(eff), n_strata = .N, consistent = uniqueN(sign(eff)) == 1),
          by = .(set, .pvn)]
  E[, size := sapply(as.character(set), function(x) length(allsets[[x]]))]
  E[, bin := szbin(size)]; E[, is_null := grepl("^__null_", as.character(set))]
  nl <- E[is_null == TRUE, .(mean0 = mean(eff), sd0 = sd(eff), n_null = .N), by = .(.pvn, bin)]
  E <- merge(E, nl, by = c(".pvn", "bin"), all.x = TRUE)
  E[, z := (eff - mean0) / sd0]; setnames(E, ".pvn", "pvn")
  nulz <- E[is_null == TRUE, .(pvn, bin, zn = z)]
  E[, p_emp := mapply(function(pv, bn, zz) {
    d <- nulz[pvn == pv & bin == bn, zn]
    (1 + sum(abs(d) >= abs(zz), na.rm = TRUE)) / (length(d) + 1) }, pvn, bin, z)]
  E[, p_norm := 2 * pnorm(-abs(z))]
  list(E = E, nulls = nl)
}

o_all <- niche_eff(L); E <- o_all$E; nulls <- o_all$nulls
if (any(is_dn(names(allsets)))) {
  o_ho <- niche_eff(L[.half == "test"])
  E <- rbind(E[!is_dn(set)], o_ho$E[is_dn(set)])
}
E[, held_out := is_dn(set)]
res <- E[is_null == FALSE][
  , .(program = as.character(set), level = pvn, n_genes = size,
      effect = round(eff, 5), z = round(z, 2), consistent,
      basis = fifelse(is_dn(set), "de novo (descriptive)", "curated"),
      p = signif(p_norm, 3), p_emp, p_norm)]
res[, q := NA_real_]
res[basis == "curated", q := signif(p.adjust(p_norm, "BH"), 3)]
res[, called := !is.na(q) & abs(z) >= Z_GATE & q <= Q_GATE & consistent]
res[, p_norm := NULL]; setorder(res, -called, p_emp, -z)

# ==== confound diagnostics =================================================
cm <- as.data.table(so@meta.data)[, .(pvn = .pvn, sec = .sec, reg = .reg,
                                      ndet = get(NDET_COL))]
conf <- cm[, .(n_cells = .N, med_detected = as.numeric(median(ndet)),
               n_sections = uniqueN(sec),
               top_region = names(sort(table(reg), decreasing = TRUE))[1],
               region_purity = round(max(table(reg)) / .N, 2)), by = pvn][order(pvn)]
conf[, region_confounded := has_region & region_purity >= 0.90]
fwrite(conf, file.path(OUT_RES, "EC_PVN_program_confounds.csv"))
print(conf)

# ==== region main effect (if region present) ==============================
if (has_region) {
  region_eff <- function(Lx) {
    gr2 <- Lx[, .(m = mean(s), n = .N), by = .(set, .stratum, .reg)][n >= MIN_CELLS]
    gr2[, sec := tstrsplit(.stratum, "|", fixed = TRUE)[[1]]]
    gr2[, eff := m - mean(m), by = .(set, sec)]
    E2 <- gr2[, .(eff = mean(eff), consistent = uniqueN(sign(eff)) == 1), by = .(set, .reg)]
    E2[, `:=`(size = sapply(as.character(set), function(x) length(allsets[[x]])),
              is_null = grepl("^__null_", as.character(set)))]
    E2[, bin := szbin(size)]
    n2 <- E2[is_null == TRUE, .(mean0 = mean(eff), sd0 = sd(eff)), by = .(.reg, bin)]
    E2 <- merge(E2, n2, by = c(".reg", "bin"), all.x = TRUE)[, z := (eff - mean0) / sd0]
    nz2 <- E2[is_null == TRUE, .(.reg, bin, zn = z)]
    E2[, p_emp := mapply(function(rg, bn, zz) {
      d <- nz2[.reg == rg & bin == bn, zn]
      (1 + sum(abs(d) >= abs(zz), na.rm = TRUE)) / (length(d) + 1) }, .reg, bin, z)]
    E2[, p_norm := 2 * pnorm(-abs(z))]; E2
  }
  E2 <- region_eff(L)
  if (any(is_dn(names(allsets))))
    E2 <- rbind(E2[!is_dn(set)], region_eff(L[.half == "test"])[is_dn(set)])
  r2 <- E2[is_null == FALSE][
    , .(program = as.character(set), level = paste0("REGION:", .reg), n_genes = size,
        effect = round(eff, 5), z = round(z, 2), consistent,
        basis = fifelse(is_dn(set), "de novo (descriptive)", "curated"),
        p = signif(p_norm, 3), p_emp, p_norm)]
  r2[, q := NA_real_]; r2[basis == "curated", q := signif(p.adjust(p_norm, "BH"), 3)]
  r2[, p_norm := NULL]
  r2[, called := !is.na(q) & abs(z) >= Z_GATE & q <= Q_GATE & consistent]
  res <- rbind(res, r2)
}

# ==== calibration ==========================================================
cal <- E[is_null == TRUE, .(n = .N, sd_z = round(sd(z, na.rm = TRUE), 2),
                            pct_abs_z_over_gate = round(100 * mean(abs(z) >= Z_GATE, na.rm = TRUE), 2),
                            pct_would_be_called = round(100 * mean(abs(z) >= Z_GATE & consistent,
                                                                   na.rm = TRUE), 2)), by = pvn]
sd0m <- nulls[, .(sd0 = mean(sd0)), by = .pvn]; setnames(sd0m, ".pvn", "pvn")
cal <- merge(cal, sd0m, by = "pvn")
cal <- merge(cal, conf[, .(pvn, med_detected, n_cells)], by = "pvn")
cal[, precision_flag := fifelse(sd0 > 1.5 * median(sd0), "LOW PRECISION", "ok")]
setnames(cal, "pvn", "level")
fwrite(cal, file.path(OUT_RES, "EC_PVN_program_calibration.csv"))
say("\n=== calibration (want sd_z ~= 1 and pct_would_be_called << 1) ===")
print(cal[, .(level, n_cells, med_detected, sd_z, pct_would_be_called,
              sd0 = round(sd0, 4), precision_flag)])
expect_pct <- 100 * 2 * pnorm(-Z_GATE)
cal[, tail_vs_normal := round(pct_abs_z_over_gate / expect_pct, 1)]
bad <- cal[sd_z > 1.4 | sd_z < 0.7 | tail_vs_normal > 3]
if (nrow(bad)) say("WARNING: normal tail not honest for ", paste(bad$level, collapse = ", "),
                   " - use p_emp and treat calls as leads.")

# reciprocal "shadow" calls (a program up in one niche is mechanically down elsewhere)
res[, shadow := FALSE]
nres <- res[!grepl("^REGION:", level)]
if (nres[called == TRUE, .N]) {
  top <- nres[called == TRUE, .(zmax = z[which.max(abs(z))]), by = program]
  res[nres[top, on = "program"][called == TRUE & sign(z) != sign(zmax) & abs(z) < abs(zmax),
                                .(program, level)], on = .(program, level), shadow := TRUE]
}
fwrite(res, file.path(OUT_RES, "EC_PVN_programs.csv"))
say("\n=== programs called (|z| >= ", Z_GATE, ", q <= ", Q_GATE, ", consistent) ===")
if (res[called == TRUE, .N]) print(res[called == TRUE][
  , .(program, level, n_genes, z, q, shadow)]) else
    say("none called at this noise floor - see EC_PVN_program_calibration.csv")

# ==== figure ===============================================================

programs_plot <- rev(c(
  'IpEC',
  'endothelial_adhesion',  "inflammatory_activation" ,
  'bbb_transport',
  'tumour_btb', #vascular remodelling
  'endmt',
  'tip_cell',
  'stalk_quiescent',
  'proliferating',
  'ifn_response',
  'ifn_antigen',
  'leukocyte_recruitment',
  'hev'
))

H <- res[!grepl("^REGION:", level) & !grepl("^denovo", program)]
H <- H |> 
  filter(program %in% programs_plot) |> 
  mutate(level = str_remove_all(level , '_'),
         level = factor(level, levels = paste0('PVN', 1:length(unique(H$level)))),
         program = factor(program, levels = programs_plot))
H <- as.data.table(H)
if (nrow(H)) {
  H[, star := fifelse(called, "*", "")]
  lim <- max(3.5, quantile(abs(H$z), 0.99, na.rm = TRUE))
  H[, zc := pmax(pmin(z, lim), -lim)]
  p1 <- ggplot(H, aes(level, program, fill = zc)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    #geom_text(aes(label = star), size = 4, vjust = 0.75, colour = "grey15") +
    scale_fill_gradient2(low = "#2166AC", 
                         mid = "white", 
                         high = "gold",
                         midpoint = 0, 
                         limits = c(-5, 5),
                         oob = scales::squish,
                         name = "z vs\nmatched\nrandom sets") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  print(p1)
  
  svglite::svglite(file.path(OUT_FIG, "pvn_heatmap_signatures.svg"), 
                   width = 2.4*1.5, height = 2.2*1.1)
  print(p1)
  dev.off()
  
#  cal2 <- copy(cal)[!grepl("^REGION", level)]
#  p2 <- ggplot(cal2, aes(reorder(level, sd0), sd0, fill = precision_flag)) +
#    geom_col(width = 0.65) +
#    scale_fill_manual(values = c("ok" = "grey65", "LOW PRECISION" = "#B23A48"), name = NULL) +
#    coord_flip() +
#    labs(x = NULL, y = "null spread of niche effect (sd0)", title = "Noise floor per niche",
#         subtitle = "larger sd0 = less sensitive; driven by cell count and detection depth") +
#    theme_bw(base_size = 9) + theme(legend.position = "top")
#  ggsave(file.path(OUT_FIG, "EC_PVN_programs.pdf"),
#         p1 / p2 + plot_layout(heights = c(2.4, 1)),
#         width = 9, height = max(6, 1.6 + 0.22 * uniqueN(H$program) + 2.2), limitsize = FALSE)
}
say("\nWrote: EC_PVN_programs.csv, EC_PVN_program_calibration.csv, ",
    "EC_PVN_program_confounds.csv; figure in ", OUT_FIG,
    " -> results in ", OUT_RES)
