################################################################################
# 03_liana_pvn_scoring.R
# Project: Perivascular invasion paper
# Author: Simon Castillo
# Last update: 19/09/2026
################################################################################
##
##


suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(spatstat.geom)
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
})

# Resolve paths independently of the current working directory.
.args_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(.args_file)) stop("Run this script with Rscript")
.script_path <- normalizePath(sub("^--file=", "", .args_file[[1]]), mustWork = TRUE)
PROJECT_ROOT <- dirname(dirname(.script_path))
.config_file <- Sys.getenv("LIANA_CONFIG", file.path(PROJECT_ROOT, "config", "config.R"))
if (!file.exists(.config_file)) stop("Configuration not found: ", .config_file)
source(.config_file, local = FALSE)

# Script-specific parameters are set in config/config.R.
SUBSETS <- MODES
CONSENSUS_FILE <- "liana_consensus.csv"
NODE_OF <- c(Endothelial="E", AC_like="T", MES_like="T",
             NPC_like="T", OPC_like="T", Macrophage_BAM="M", Microglia_core="M")
OUT_DIR <- PVN_LR_OUT
N_PERM <- PVN_N_PERM
SEED <- SEED_PVN
if (length(QUPATH) != 2L) stop("PVN reproducibility table currently requires exactly two images")
if (!identical(names(QUPATH), names(EXPR))) stop("QUPATH/EXPR image keys differ")
if (!file.exists(PVN_ASSIGN)) stop("Missing PVN assignment CSV: ", PVN_ASSIGN)
if (!dir.exists(CONSENSUS_DIR)) stop("Missing consensus directory: ", CONSENSUS_DIR)
dir.create(dirname(GENE_AXIS), recursive = TRUE, showWarnings = FALSE)

AXIS_RULES <- list(
  vascular_activation = c("^Vegf", "^Kdr$", "^Flt[14]$", "^Angpt", "^Tek$", "^Tie1$",
                          "^Esm1$", "^Apln", "^Dll[34]$", "^Notch[1-4]$", "^Jag[12]$",
                          "^Efnb2$", "^Ephb4$", "^Eng$", "^Acvrl1$", "^Cdh5$",
                          "^Pecam1$", "^Ptprb$", "^Robo4$", "^Unc5b$", "^S1pr[1-3]$",
                          "^Edn[123]$", "^Ednr[ab]$", "^Nos3$", "^Sel[ep]$",
                          "^Vcam1$", "^Icam[12]$", "^Pdgfb$", "^Pdgfrb$", "^Nrp[12]$",
                          "^Cldn5$", "^Slc2a1$", "^Mfsd2a$", "^Plvap$", "^Ackr1$",
                          "^Cav1$", "^Emcn$", "^Podxl$", "^Sox1[78]$", "^Esam$",
                          "^F11r$", "^Jam[23]$", "^Mcam$", "^Cd93$", "^Pigf$"),
  ecm_adhesion        = c("^Col[0-9]", "^Lam[abc][0-9]", "^Fn1$", "^Agrn$", "^Nid[12]$",
                          "^Hspg2$", "^Tnc$", "^Tnr$", "^Thbs[1-4]$", "^Spp1$", "^Vtn$",
                          "^Sdc[1-4]$", "^Itga", "^Itgb", "^Dag1$", "^Cd44$", "^Mmp",
                          "^Timp", "^Bgn$", "^Dcn$", "^Lum$", "^Postn$", "^Fbn[12]$",
                          "^Vcan$", "^Bcan$", "^Ncan$", "^Emilin", "^Sparc$", "^Fbln",
                          "^Ecm1$", "^Lgals[13]$", "^Anxa2$", "^Plaur$", "^Plau$",
                          "^Plat$", "^Ddr1$", "^Vcl$", "^Tln1$", "^Ltbp[13]$", "^Has2$",
                          "^Edil3$", "^Mfge8$", "^Cspg4$", "^Serpine[12]$"),
  immune_antigen      = c("^H2-", "^B2m$", "^Cd74$", "^Cd4$", "^Cd8[ab]", "^Ciita$",
                          "^Tap[12]$", "^Ccl[0-9]", "^Ccr[0-9]", "^Cxcl[0-9]",
                          "^Cxcr[0-9]", "^Cx3cl1$", "^Il[0-9]", "^Il1r", "^Tnfsf",
                          "^Tnfrsf", "^C1q", "^C3$", "^C3ar1$", "^C4b$", "^Cfh$",
                          "^Serping1$", "^Tlr[124]$", "^Cd180$", "^Ly86$", "^Cd80$",
                          "^Cd86$", "^Ctla4$", "^Pdcd1", "^Cd274$", "^Trem2$",
                          "^Tyrobp$", "^Csf[123]", "^Csf1r$", "^Csf2ra$", "^Csf3r$",
                          "^Itga[mlx]$", "^Itgb2$", "^Ptprc$", "^Lgals9$", "^Havcr2$",
                          "^Sirpa$", "^Cd47$", "^Cd68$", "^Mrc1$", "^Fcgr", "^Clec12a$",
                          "^Cd300", "^Selplg$", "^P2ry12$", "^Gpnmb$", "^Cd36$",
                          "^Mertk$", "^Gas6$", "^Axl$", "^Tyro3$", "^Pros1$", "^Mif$",
                          "^Hsp90", "^Calr$", "^S100a"),
  trophic_neural      = c("^Ptn$", "^Ptprz1$", "^Mdk$", "^Bdnf$", "^Ntrk[1-3]$", "^Ngf",
                          "^Gdnf$", "^Gfra1$", "^Nrg[1-4]$", "^Erbb[2-4]$", "^Sema[0-9]",
                          "^Plxn", "^Nrxn[1-3]$", "^Nlgn[1-3]$", "^Cntn", "^Ncam[12]$",
                          "^L1cam$", "^Chl1$", "^Nrcam$", "^Nfasc$", "^Dcc$", "^Neo1$",
                          "^Ntn[1-4]$", "^Robo[1-3]$", "^Slit[1-3]$", "^Slitrk", "^Efna",
                          "^Epha[0-9]", "^Unc5[ac]$", "^Rgm[ab]$", "^Rtn4", "^Omg$",
                          "^Lingo1$", "^Mag$", "^Gria[1-4]$", "^Grin[12]", "^Gabbr[12]$",
                          "^Atp1a[1-3]$", "^Igf[12]", "^Igf1r$", "^Insr$", "^Fgf[0-9]",
                          "^Fgfr[1-4]$", "^Egfr$", "^Hbegf$", "^Areg$", "^Wnt[0-9]",
                          "^Fzd[0-9]", "^Lrp[4568]$", "^App$", "^Aplp[12]$", "^Bmp[0-9]",
                          "^Bmpr", "^Acvr", "^Tgfb", "^Tgfbr", "^Gpc[1-6]$", "^Nptx",
                          "^Cadm[1-4]$", "^Gpr37", "^Dlk1$", "^Lifr$"))

set.seed(SEED)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

find_col <- function(nms, patterns, what, required = TRUE) {
  for (p in patterns) {
    hit <- grep(p, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  if (required)
    stop(sprintf("Could not find the %s column. Available: %s", what,
                 paste(nms, collapse = ", ")), call. = FALSE)
  NA_character_
}

node_of <- function(ct) {
  n <- unname(NODE_OF[ct]); n[is.na(n)] <- ct[is.na(n)]; n
}

# ---- 1. PVN assignments ---------------------------------------------------
assign_dt <- fread(PVN_ASSIGN)
stopifnot(all(c("Image", "cell_id", "PVN") %in% names(assign_dt)))
if (REPLICATING_ONLY && "replicating" %in% names(assign_dt)) {
  n0 <- nrow(assign_dt)
  assign_dt <- assign_dt[replicating == TRUE]
  cat(sprintf("Anchors: %d of %d in replicating PVNs\n", nrow(assign_dt), n0))
}
assign_dt[, PVN := if (length(PVN_LABELS)) coalesce(PVN_LABELS[PVN], PVN) else PVN]
cat("\n== ANCHORS PER PVN ==\n")
print(dcast(assign_dt[, .N, by = .(Image, PVN)], PVN ~ Image, value.var = "N", fill = 0))

# ---- 2. consensus LR pairs from the three subset files -------------------
read_consensus <- function(sub) {
  f <- file.path(CONSENSUS_DIR, sub, CONSENSUS_FILE)
  if (!file.exists(f)) stop("Missing consensus file: ", f)
  x <- fread(f)
  lc <- find_col(names(x), c("^ligand\\.complex$", "^ligand$", "ligand"),   "ligand")
  rc <- find_col(names(x), c("^receptor\\.complex$", "^receptor$", "receptor"), "receptor")
  sc <- find_col(names(x), c("^source$", "^sender$", "source"), "source")
  tc <- find_col(names(x), c("^target$", "^receiver$", "target"), "target")
  qc <- find_col(names(x), c("^q_fisher$", "qvalue", "q_value", "^padj$", "fdr",
                             "^q$", "aggregate_rank"), "q (significance)", required = FALSE)
  mc <- find_col(names(x), c("^mean_LRscore$", "LRscore", "magnitude", "consensus_score"),
                 "magnitude", required = FALSE)
  out <- data.table(subset = sub,
                    source = x[[sc]], target = x[[tc]],
                    ligand.complex = x[[lc]], receptor.complex = x[[rc]],
                    q_use = if (!is.na(qc)) as.numeric(x[[qc]]) else NA_real_,
                    mean_LRscore = if (!is.na(mc)) as.numeric(x[[mc]]) else NA_real_)
  cat(sprintf("  %s: %d pairs (%s: %s)\n", sub, nrow(out),
              if (is.na(qc)) "no q col" else qc,
              if (is.na(mc)) "no magnitude col" else mc))
  out
}
cat("\n== reading split consensus ==\n")
cons <- rbindlist(lapply(SUBSETS, read_consensus))
stopifnot(nrow(cons) > 0)

if (any(!is.na(cons$q_use))) {
  n0 <- nrow(cons); cons <- cons[is.na(q_use) | q_use <= Q_MAX]
  cat(sprintf("Kept %d of %d pairs at q<=%.3g\n", nrow(cons), n0, Q_MAX))
} else {
  message("No significance column found in any consensus file - Q_MAX ignored, all pairs kept.")
}

cons[, src_node := node_of(source)][, tgt_node := node_of(target)]
if (EC_ONLY) {
  n0 <- nrow(cons); cons <- cons[source == EC_CLASS | target == EC_CLASS]
  cat(sprintf("EC_ONLY: kept %d of %d EC-involving pairs\n", nrow(cons), n0))
}
cons[, node_edge := paste(src_node, tgt_node, sep = "->")]
cons[, direction := fifelse(source == EC_CLASS & target == EC_CLASS, "autocrine_EC",
                            fifelse(source == EC_CLASS, "from_EC",
                                    fifelse(target == EC_CLASS, "to_EC", node_edge)))]
cons[, partner := fifelse(source == EC_CLASS, target, source)]
cons[, pair_id := paste(source, target, ligand.complex, receptor.complex, sep = "|")]
# a pair reported in two subsets (e.g. EC autocrine in E_T and M_E) -> keep best q
setorder(cons, q_use, na.last = TRUE)
cons <- unique(cons, by = "pair_id")
cat(sprintf("\n%d unique consensus pairs after de-duplication\n", nrow(cons)))
cat("== pairs per node edge ==\n")
print(cons[, .N, by = node_edge][order(-N)])

# ---- 3. biological axes ---------------------------------------------------
# seed the per-gene table over the consensus gene universe if it does not exist
if (!file.exists(GENE_AXIS)) {
  gu <- sort(unique(unlist(strsplit(
    c(cons$ligand.complex, cons$receptor.complex), "[_:]"))))
  axis_of_gene <- function(s) {
    hit <- names(AXIS_RULES)[vapply(AXIS_RULES, function(ps)
      any(vapply(ps, function(p) grepl(p, s), logical(1))), logical(1))]
    if (!length(hit)) "unassigned" else paste(sort(hit), collapse = "+")
  }
  ga <- data.table(gene = gu, axis = vapply(gu, axis_of_gene, character(1)))
  ga[, needs_review := grepl("\\+|unassigned", axis)]
  fwrite(ga, GENE_AXIS)
  cat(sprintf("\nSEEDED %s with %d genes (%d need review). CURATE it, then re-run.\n",
              GENE_AXIS, nrow(ga), sum(ga$needs_review)))
}
gene_axis <- as.data.table(read_csv(GENE_AXIS, show_col_types = FALSE))
stopifnot(all(c("gene", "axis") %in% names(gene_axis)))
axis_lookup <- setNames(gene_axis$axis, gene_axis$gene)

axis_of <- function(sym) {
  parts <- strsplit(sym, "[_:]")[[1]]
  hits  <- unique(unname(axis_lookup[parts])); hits <- hits[!is.na(hits)]
  if (!length(hits)) "unassigned" else paste(sort(hits), collapse = "+")
}

cover <- setdiff(unlist(strsplit(unique(c(cons$ligand.complex, cons$receptor.complex)),
                                 "[_:]")), names(axis_lookup))
if (length(cover))
  message(sprintf("%d consensus symbols are not in %s and will be 'unassigned': %s",
                  length(cover), GENE_AXIS, paste(head(cover, 20), collapse = ", ")))

## cached AXIS_FILE wins, but only while newer than the tables it was built from
if (file.exists(AXIS_FILE)) {
  srcs <- c(GENE_AXIS, COMBO_AXIS)[file.exists(c(GENE_AXIS, COMBO_AXIS))]
  if (length(srcs) && file.mtime(AXIS_FILE) < max(file.mtime(srcs))) {
    cat(sprintf("%s is older than %s - regenerating from source tables.\n",
                AXIS_FILE, paste(basename(srcs), collapse = "/")))
    file.rename(AXIS_FILE, paste0(AXIS_FILE, ".superseded"))
  }
}
if (file.exists(AXIS_FILE)) {
  axis_map <- as.data.table(read_csv(AXIS_FILE, show_col_types = FALSE))
  miss <- setdiff(cons$pair_id, axis_map$pair_id)
  if (length(miss))
    stop(sprintf("%d consensus pairs are missing from %s. Delete it to regenerate.",
                 length(miss), AXIS_FILE), call. = FALSE)
  if (!"sign" %in% names(axis_map)) axis_map[, sign := NA_character_]
  cat(sprintf("Using curated axis assignment (%s)\n", AXIS_FILE))
} else {
  axis_map <- copy(cons)[, .(pair_id, source, target, ligand.complex,
                             receptor.complex, direction, partner)]
  axis_map[, axis_ligand   := vapply(ligand.complex,   axis_of, character(1))]
  axis_map[, axis_receptor := vapply(receptor.complex, axis_of, character(1))]
  axis_map[, axis := fifelse(axis_receptor == "unassigned", axis_ligand, axis_receptor)]
  
  ## ---- position-aware exclusion --------------------------------------------
  ## The RECEIVER decides whether a pair is scoreable. An 'exclude' (non-surface)
  ## receptor cannot receive signal -> the pair is dropped. An 'exclude' LIGAND is
  ## not fatal: many excluded genes are moonlighting / DAMP ligands (e.g. Calr,
  ## Hsp90aa1/b1, Vim, Snca, S100a1) that are secreted or surface-exposed by
  ## stressed/dying cells, so such a pair is KEPT and scored on the receptor's
  ## axis. Set DROP_IF_LIGAND_EXCLUDED <- TRUE (anywhere above) to restore the old
  ## both-sides behaviour.
  if (DROP_IF_LIGAND_EXCLUDED) {
    axis_map[axis_ligand == "exclude" | axis_receptor == "exclude", axis := "exclude"]
  } else {
    axis_map[axis_receptor == "exclude", axis := "exclude"]
    n_resc <- axis_map[axis_ligand == "exclude" &
                         !axis_receptor %in% c("exclude", "unassigned"), .N]
    cat(sprintf(paste0("Position-aware exclusion: kept %d pair(s) with an 'exclude' ",
                       "ligand but a real receptor axis (scored on the receptor axis; ",
                       "flagged needs_review).\n"), n_resc))
  }
  ## --------------------------------------------------------------------------
  axis_map[, needs_review := axis_ligand != axis_receptor |
             grepl("\\+|unassigned", paste(axis_ligand, axis_receptor))]
  axis_map[, lr_combo := paste(ligand.complex, receptor.complex, sep = "->")]
  axis_map[, sign := NA_character_]
  if (file.exists(COMBO_AXIS)) {
    cb  <- as.data.table(read_csv(COMBO_AXIS, show_col_types = FALSE))
    hit <- match(axis_map$lr_combo, cb$lr_combo); ok <- !is.na(hit)
    axis_map[ok, `:=`(axis = cb$axis[hit[ok]], sign = cb$sign[hit[ok]],
                      needs_review = FALSE)]
    cat(sprintf("Combo curation applied to %d of %d pairs (%s)\n",
                sum(ok), nrow(axis_map), COMBO_AXIS))
  }
  fwrite(axis_map, AXIS_FILE)
  cat(sprintf("\nWrote proposed axis assignment: %s\n", AXIS_FILE))
  cat(sprintf("  %d of %d pairs flagged needs_review - CURATE, then re-run.\n",
              sum(axis_map$needs_review), nrow(axis_map)))
  print(axis_map[, .N, by = .(axis, needs_review)][order(axis)])
}

lr <- merge(cons[, .(pair_id, source, target, ligand.complex, receptor.complex,
                     direction, partner, node_edge, mean_LRscore)],
            axis_map[, .(pair_id, axis, sign)], by = "pair_id")[!axis %in% c("unassigned", "exclude")]
lr[, lr_combo := paste(ligand.complex, receptor.complex, sep = "->")]
cat(sprintf("Scoring %d of %d consensus pairs (%d unassigned, %d excluded)\n",
            nrow(lr), nrow(cons), sum(axis_map$axis == "unassigned"),
            sum(axis_map$axis == "exclude")))

# ---- 4. per-image scoring -------------------------------------------------
score_image <- function(img, keep_pairs = NULL) {
  d <- fread(QUPATH[[img]])
  setnames(d, grep("^Centroid X", names(d), value = TRUE), "cx")
  setnames(d, grep("^Centroid Y", names(d), value = TRUE), "cy")
  id_q <- find_col(names(d), c("^cell_id$", "^cell$", "^Object ID$"), "cell id (QuPath)")
  d <- d[!Classification %in% c("Ependymal", "not_assigned", "multiple")]
  d[, cid := as.character(get(id_q))]
  
  obj  <- readRDS(EXPR[[img]])
  mat  <- obj$matrix
  meta <- as.data.table(obj$metadata)
  id_m <- find_col(names(meta), c("^cell_id$", "^cell$", "^barcode$"), "cell id (expr)")
  meta[, cid := as.character(get(id_m))]
  meta <- meta[match(colnames(mat), cid)]
  stopifnot(identical(meta$cid, colnames(mat)))
  
  common <- intersect(d$cid, meta$cid)
  cat(sprintf("\n[%s] %d QuPath cells, %d expression cells, %d joined by cell id (%.1f%% of QuPath)\n",
              img, nrow(d), ncol(mat), length(common), 100 * length(common) / nrow(d)))
  if (length(common) < 0.5 * nrow(d))
    stop("Fewer than half the QuPath cells matched the expression matrix - the ",
         "cell id columns are probably not the same key.", call. = FALSE)
  d <- d[cid %in% common]
  d[, label := meta[[LABEL_COL]][match(cid, meta$cid)]]
  d[, col_i := match(cid, colnames(mat))]
  d[, gid := .I]
  
  anchors <- copy(assign_dt[Image == img])
  anchors[, cid_chr := as.character(cell_id)]
  d[, is_anchor := cid %in% anchors$cid_chr]
  d[, PVN := anchors$PVN[match(cid, anchors$cid_chr)]]
  if ("region" %in% names(anchors))
    d[, region := anchors$region[match(cid, anchors$cid_chr)]]
  n_anchor_lost <- sum(!anchors$cid_chr %in% d$cid)
  if (n_anchor_lost)
    message(sprintf("[%s] %d of %d PVN anchors have no expression profile - dropped",
                    img, n_anchor_lost, nrow(anchors)))
  anc <- d[is_anchor == TRUE & !is.na(PVN)]
  stopifnot(nrow(anc) > 0)
  
  pp <- ppp(d$cx, d$cy, window = owin(range(d$cx), range(d$cy)), marks = d$label)
  cp <- closepairs(pp, r = R_NB, twice = FALSE)
  I <- d$gid[cp$i]; J <- d$gid[cp$j]
  a_set <- anc$gid
  keep <- I %in% a_set | J %in% a_set
  I <- I[keep]; J <- J[keep]
  edges <- data.table(anchor = fifelse(I %in% a_set, I, J),
                      nb     = fifelse(I %in% a_set, J, I))
  both <- edges[anchor %in% a_set & nb %in% a_set]
  edges <- rbind(edges, both[, .(anchor = nb, nb = anchor)])
  edges <- edges[anchor %in% a_set]
  edges[, nb_label := d$label[nb]]
  edges[, nb_col := d$col_i[nb]]
  cat(sprintf("[%s] %d anchors, %d anchor-neighbour edges, median %d neighbours per anchor\n",
              img, nrow(anc), nrow(edges), median(edges[, .N, by = anchor]$N)))
  
  genes <- unique(c(lr$ligand.complex, lr$receptor.complex))
  resolve <- function(sym) {
    if (sym %in% rownames(mat)) return(sym)
    parts <- strsplit(sym, "[_:]")[[1]]
    if (length(parts) > 1 && all(parts %in% rownames(mat))) return(parts)
    NA_character_
  }
  resolved <- setNames(lapply(genes, resolve), genes)
  absent <- names(resolved)[vapply(resolved, function(z) all(is.na(z)), logical(1))]
  if (length(absent))
    message(sprintf("[%s] %d LR symbols absent from the matrix: %s",
                    img, length(absent), paste(head(absent, 15), collapse = ", ")))
  used <- setdiff(names(resolved), absent)
  sub  <- mat[unique(unlist(resolved[used])), , drop = FALSE]
  if (SCALE == "geomean_linear") {
    sub@x <- expm1(sub@x)
  }
  ex <- do.call(rbind, lapply(used, function(sym) {
    g <- resolved[[sym]]
    if (length(g) == 1) sub[g, , drop = FALSE] else
      Matrix(apply(as.matrix(sub[g, , drop = FALSE]), 2, COMPLEX_FN), nrow = 1,
             sparse = TRUE)
  }))
  rownames(ex) <- used
  
  anc_idx <- setNames(seq_len(nrow(anc)), anc$gid)
  classes <- sort(unique(edges$nb_label))
  S <- lapply(classes, function(cl) {
    e <- edges[nb_label == cl]
    M <- sparseMatrix(i = e$nb_col, j = anc_idx[as.character(e$anchor)], x = 1,
                      dims = c(ncol(ex), nrow(anc)))
    list(sum = ex %*% M, n = Matrix::colSums(M))
  })
  names(S) <- classes
  A <- ex[, anc$col_i, drop = FALSE]
  
  pvn_f <- factor(anc$PVN)
  G0 <- sparseMatrix(i = seq_len(nrow(anc)), j = as.integer(pvn_f), x = 1,
                     dims = c(nrow(anc), nlevels(pvn_f)),
                     dimnames = list(NULL, levels(pvn_f)))
  
  pooled <- function(G) {
    n_anc <- Matrix::colSums(G)
    ec <- as.matrix(A %*% G) %*% diag(1 / n_anc, nrow = length(n_anc))
    dimnames(ec) <- list(rownames(ex), colnames(G))
    nb <- lapply(S, function(s) {
      cnt <- as.numeric(s$n %*% G)
      m <- as.matrix(s$sum %*% G)
      m <- sweep(m, 2, pmax(cnt, 1), "/")
      m[, cnt < MIN_CELLS_PVN] <- 0
      dimnames(m) <- list(rownames(ex), colnames(G))
      list(mean = m, n = cnt)
    })
    list(ec = ec, nb = nb, n_anchor = n_anc)
  }
  
  lr_i <- lr[ligand.complex %in% used & receptor.complex %in% used]
  if (!is.null(keep_pairs)) lr_i <- lr_i[pair_id %in% keep_pairs]
  if (!nrow(lr_i)) stop(sprintf("[%s] no pairs left after common-support filter", img))
  score_from <- function(P) {
    vapply(seq_len(nrow(lr_i)), function(k) {
      r <- lr_i[k]
      lig <- if (r$source == EC_CLASS) P$ec[r$ligand.complex, ] else {
        cl <- P$nb[[r$source]]; if (is.null(cl)) rep(0, ncol(P$ec)) else
          cl$mean[r$ligand.complex, ] }
      rec <- if (r$target == EC_CLASS) P$ec[r$receptor.complex, ] else {
        cl <- P$nb[[r$target]]; if (is.null(cl)) rep(0, ncol(P$ec)) else
          cl$mean[r$receptor.complex, ] }
      if (SCALE == "linear") lig * rec else sqrt(lig * rec)
    }, numeric(ncol(P$ec)))
  }
  
  P0  <- pooled(G0)
  obs <- score_from(P0)
  dimnames(obs) <- list(colnames(G0), lr_i$pair_id)
  
  ## which pair x PVN cells were zero-filled because a NON-EC partner class was
  ## rare/absent among that PVN's neighbours. Generalised from the original:
  ## a myeloid<->tumour pair needs BOTH classes present, so it is flagged if
  ## EITHER falls below the floor.
  low <- vapply(seq_len(nrow(lr_i)), function(k) {
    need <- setdiff(c(lr_i$source[k], lr_i$target[k]), EC_CLASS)
    if (!length(need)) return(rep(FALSE, ncol(G0)))   # pure EC autocrine
    Reduce(`|`, lapply(need, function(cl) {
      nb <- P0$nb[[cl]]
      if (is.null(nb)) rep(TRUE, ncol(G0)) else nb$n < MIN_CELLS_PVN
    }))
  }, logical(ncol(G0)))
  dimnames(low) <- dimnames(obs)
  
  axis_of_pair <- lr_i$axis
  axis_pairs   <- split(seq_len(nrow(lr_i)), axis_of_pair)
  axis_I <- function(M)
    t(vapply(axis_pairs, function(j) rowSums(M[, j, drop = FALSE]), numeric(nrow(M))))
  
  rel_of <- function(I) {
    med <- apply(I, 1, median, na.rm = TRUE)
    eps <- pmax(med, .Machine$double.eps) * 1e-6
    log2(sweep(I + eps, 1, med + eps, "/"))
  }
  obs_I   <- axis_I(obs)
  obs_rel <- rel_of(obs_I)
  
  zsc <- function(M) {
    mu  <- colMeans(M, na.rm = TRUE)
    sdv <- apply(M, 2, sd, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv == 0] <- NA_real_
    z <- sweep(sweep(M, 2, mu, "-"), 2, sdv, "/")
    z[!is.finite(z)] <- 0
    z
  }
  combo_idx    <- split(seq_len(nrow(lr_i)), lr_i$lr_combo)
  combo_axis   <- vapply(combo_idx, function(j) lr_i$axis[j[1]], character(1))
  axis_combos  <- split(names(combo_idx), combo_axis)
  axis_score <- function(L) {
    cm <- vapply(combo_idx, function(j)
      rowMeans(L[, j, drop = FALSE], na.rm = TRUE), numeric(nrow(L)))
    if (!is.matrix(cm)) cm <- matrix(cm, nrow = nrow(L),
                                     dimnames = list(NULL, names(combo_idx)))
    t(vapply(axis_combos, function(cs)
      rowMeans(cm[, cs, drop = FALSE], na.rm = TRUE), numeric(nrow(L))))
  }
  obs_bal <- axis_score(zsc(obs))[rownames(obs_rel), , drop = FALSE]
  
  perm <- array(NA_real_, c(nrow(obs_rel), ncol(obs_rel), N_PERM),
                dimnames = list(rownames(obs_rel), colnames(obs_rel), NULL))
  permb <- perm
  strata <- if (PERM_WITHIN_REGION && "region" %in% names(anc))
    as.integer(factor(anc$region)) else rep(1L, nrow(anc))
  if (PERM_WITHIN_REGION) cat(sprintf("[%s] permuting within %d region strata\n",
                                      img, length(unique(strata))))
  for (b in seq_len(N_PERM)) {
    ii <- ave(seq_len(nrow(anc)), strata, FUN = function(z) sample(z))
    Gb <- G0[ii, , drop = FALSE]
    sb <- score_from(pooled(Gb))
    perm[, , b]  <- rel_of(axis_I(sb))
    permb[, , b] <- axis_score(zsc(sb))[rownames(obs_rel), , drop = FALSE]
  }
  pval <- function(A, o) {
    hit <- apply(abs(A) >= abs(as.vector(o)), c(1, 2), sum, na.rm = TRUE)
    den <- apply(!is.na(A), c(1, 2), sum)
    if (min(den) < 0.8 * N_PERM)
      message(sprintf("[%s] some axes lost permutations to the cell floor (min %d of %d usable)",
                      img, min(den), N_PERM))
    (1 + hit) / (1 + den)
  }
  p_two <- pval(perm,  obs_rel)
  p_bal <- pval(permb, obs_bal)
  
  axis_dt <- as.data.table(as.table(obs_rel))
  setnames(axis_dt, c("axis", "PVN", "axis_log2fc"))
  axis_dt[, I_raw   := as.vector(as.table(obs_I))]
  axis_dt[, axis_z_balanced := as.vector(as.table(obs_bal))]
  axis_dt[, p_perm := as.vector(as.table(p_two))]
  axis_dt[, p_perm_balanced := as.vector(as.table(p_bal))]
  axis_dt[, null_sd := as.vector(apply(perm, c(1, 2), sd, na.rm = TRUE))]
  axis_dt[, n_pairs  := as.integer(table(axis_of_pair)[axis])]
  axis_dt[, n_combos := lengths(axis_combos)[axis]]
  axis_dt[, thin_axis := n_pairs < MIN_PAIRS_PER_AXIS | n_combos < MIN_COMBOS_PER_AXIS]
  zf <- vapply(axis_pairs, function(j) rowSums(low[, j, drop = FALSE]), numeric(nrow(low)))
  axis_dt[, n_zero_filled := as.vector(as.table(t(zf)))]
  axis_dt[, q_perm := p.adjust(p_perm, "BH")]
  axis_dt[, Image := img]
  
  pair_dt <- as.data.table(as.table(obs))
  setnames(pair_dt, c("PVN", "pair_id", "score"))
  pair_dt[, z_across_pvn := as.vector(as.table(zsc(obs)))]
  pair_dt[, zero_filled := as.vector(as.table(low))]
  pair_dt[, Image := img]
  pair_dt <- merge(pair_dt, lr_i[, .(pair_id, source, target, ligand.complex,
                                     receptor.complex, direction, partner, node_edge, axis)],
                   by = "pair_id")
  
  n_dt <- rbindlist(lapply(names(P0$nb), function(cl)
    data.table(Image = img, PVN = colnames(G0), partner = cl,
               n_cells = P0$nb[[cl]]$n)))
  n_dt <- rbind(n_dt, data.table(Image = img, PVN = colnames(G0),
                                 partner = EC_CLASS, n_cells = P0$n_anchor))
  
  list(axis = axis_dt, pair = pair_dt, n = n_dt)
}

res  <- lapply(names(QUPATH), score_image)
names(res) <- names(QUPATH)
axis_sc <- rbindlist(lapply(res, `[[`, "axis"))
pair_sc <- rbindlist(lapply(res, `[[`, "pair"))
cells_n <- rbindlist(lapply(res, `[[`, "n"))

pair_sc[, pct_of_axis := if (sum(score) > 0) 100 * score / sum(score) else 0 * score,
        by = .(Image, PVN, axis)]

cat(sprintf("\nEvery PVN scored for every axis: %d PVN x axis x image rows, %d with no NA\n",
            nrow(axis_sc), sum(!is.na(axis_sc$axis_log2fc))))
cov_dt <- pair_sc[, .(pairs = .N, 
                      zero_filled = sum(zero_filled),
                      pct_zero = round(100 * mean(zero_filled), 1)),
                  by = .(Image, PVN, axis)]
fwrite(cov_dt, file.path(OUT_DIR, "pvn_axis_coverage.csv"))
cat("axis x PVN cells where >50% of pairs were zero-filled (partner class scarce):\n")
print(cov_dt[pct_zero > 50][order(-pct_zero)], nrows = 15)
top_contrib <- pair_sc[, .(pct = mean(pct_of_axis)), by = .(axis, ligand.complex,
                                                            receptor.complex)][
                                                              order(axis, -pct)][, head(.SD, 3), by = axis]
cat("\n== top contributors to each axis sum (mean % of I across PVNs) ==\n")
print(top_contrib)

fwrite(axis_sc, file.path(OUT_DIR, "pvn_axis_scores.csv"))
fwrite(dcast(axis_sc, Image + PVN ~ axis, value.var = "I_raw"),
       file.path(OUT_DIR, "pvn_axis_I_raw.csv"))
fwrite(pair_sc, file.path(OUT_DIR, "pvn_lr_pair_scores.csv"))
fwrite(cells_n, file.path(OUT_DIR, "pvn_partner_cell_counts.csv"))

# ---- 5. reproducibility across the two images ----------------------------
imgs <- names(QUPATH)
rep_dt <- dcast(axis_sc, axis + PVN ~ Image,
                value.var = c("axis_log2fc", "p_perm", "q_perm"))
lf1 <- paste0("axis_log2fc_", imgs[1]); lf2 <- paste0("axis_log2fc_", imgs[2])
q1  <- paste0("q_perm_", imgs[1]);      q2  <- paste0("q_perm_", imgs[2])
rep_dt[, same_sign := sign(get(lf1)) == sign(get(lf2))]
rep_dt[, sig_both  := get(q1) <= 0.05 & get(q2) <= 0.05]
rep_dt[, reproducible := same_sign & sig_both]
rep_dt[is.na(reproducible), reproducible := FALSE]
rep_dt[, mean_log2fc := rowMeans(cbind(get(lf1), get(lf2)))]
rep_dt[, fold_spread := pmax(abs(get(lf1)), abs(get(lf2))) /
         pmax(pmin(abs(get(lf1)), abs(get(lf2))), 1e-6)]
rep_dt[, abs_log2fc := abs(mean_log2fc)]
setorder(rep_dt, -reproducible, -abs_log2fc)
fwrite(rep_dt, file.path(OUT_DIR, "pvn_axis_reproducible.csv"))

cat(sprintf("\n%d of %d PVN x axis combinations are reproducible (same sign, q<=0.05 in both)\n",
            sum(rep_dt$reproducible), nrow(rep_dt)))
print(head(rep_dt[reproducible == TRUE, .(PVN, axis, mean_log2fc, fold_spread)], 15))


