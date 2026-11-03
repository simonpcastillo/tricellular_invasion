# =============================================================================
# PVN COMPOSITION HEATMAP  (niche x cell-class mean composition)
# =============================================================================
# The right display when niches overlap in a UMAP: it shows what each niche is
# MADE OF, and which classes distinguish it, directly. Two panels:
#   (a) mean proportion per niche (absolute composition, sequential fill)
#   (b) column z-score across niches (enrichment/depletion, diverging @ 0)
# Values are printed in every cell (< ~200 cells => readable, figure-style S6.5).
#
#   prop : proportion matrix (rows = ECs, cols = cell classes), sums ~1 per row
#   cnt  : data.table aligned to prop rows, with $PVN
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

pvn_composition_heatmap <- function(prop, cnt, view = c("both", "prop", "z"),
                                    cluster_rows = TRUE, digits = 2) {
  view <- match.arg(view)
  cls  <- colnames(prop)
  M <- as.data.table(prop); M[, PVN := cnt$PVN]
  mean_comp <- M[, lapply(.SD, mean), by = PVN, .SDcols = cls]
  nvec <- M[, .N, by = PVN][match(mean_comp$PVN, PVN), N]
  mat  <- as.matrix(mean_comp[, ..cls]); rownames(mat) <- mean_comp$PVN

  # order rows by similarity, cols by overall abundance
  if (cluster_rows && nrow(mat) > 2)
    mat <- mat[hclust(dist(mat))$order, , drop = FALSE]
  mat <- mat[, order(-colMeans(mat)), drop = FALSE]
  z   <- scale(mat)                                   # column z-score across niches
  z[!is.finite(z)] <- 0
  row_lv <- rownames(mat); col_lv <- colnames(mat)

  long <- function(m, val) {
    d <- as.data.table(as.table(m)); setnames(d, c("PVN", "class", val))
    d[, PVN   := factor(PVN, row_lv)]
    d[, class := factor(class, col_lv)]; d[]
  }
  dp <- long(mat, "prop"); dz <- long(z, "z")
  labdt <- data.table(PVN = factor(row_lv, row_lv),
                      lab = paste0(row_lv, "  (n=", nvec[match(row_lv, mean_comp$PVN)], ")"))
  ylab <- setNames(labdt$lab, as.character(labdt$PVN))

  hp <- ggplot(dp |> mutate(
    PVN = factor(str_remove_all(PVN, '_'),
                 levels = paste0("PVN", 1:length(unique(dp$PVN)))),
    proportion_bin = cut(
      prop,
      breaks = c(-Inf, 0.05, 0.1, 0.2, 0.3, Inf),
      labels = c(
        "< 0.05",
        "0.05–0.10",
        "0.10–0.20",
        "0.20–0.30",
        "> 0.30"
      ),
      include.lowest = TRUE
    )
  ), aes(class, PVN, fill = proportion_bin)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    #geom_text(aes(label = formatC(prop, format = "f", digits = digits)), size = 2.7) +
    scale_fill_brewer(
      palette = "RdPu",
      direction = 1,
      drop = FALSE) +
    scale_y_discrete(labels = ylab) +
    labs(title = "", x = NULL, y = NULL) +
    theme_minimal(base_size = 6) +
    guides(fill = 'none')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(), plot.title = element_text(face = "bold"))

  lim <- max(abs(dz$z), na.rm = TRUE)
  hz <- ggplot(dz|> mutate(
    PVN = factor(str_remove_all(PVN, '_'),
                 levels = paste0("PVN", 1:length(unique(dz$PVN)))))
    , aes(class, PVN, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    #geom_text(aes(label = formatC(z, format = "f", digits = 1)), size = 2.7) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-lim, lim), name = "z across\nniches") +
    scale_y_discrete(labels = ylab) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 6) +
    #guides(fill = 'none')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(), plot.title = element_text(face = "bold"))

  if (view == "prop") return(hp)
  if (view == "z")    return(hz)
  hp + hz
}
