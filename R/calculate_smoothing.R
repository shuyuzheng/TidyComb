# TidyComb
# Functions for calculating drug combination response or synergy surface
# Copyright Shuyu Zheng
#
# Functions on this page:
# smoothing
smoothing <-  function (scores.mat, len = 2) {

  options(scipen = 999)
  nr <- nrow(scores.mat)
  nc <- ncol(scores.mat)

  # missing value imputation
  while (sum(is.na(scores.mat))) {
    scores.mat <- ImputeNear(scores.mat)
  }
  ext.row.len <- (nr - 1) * (len + 2) - (nr - 2)
  ext.col.len <- (nc - 1) * (len + 2) - (nc - 2)

  extended.row.idx <- seq(1, nr, length = ext.row.len)
  extended.col.idx <- seq(1, nc, length = ext.col.len)

  krig.coord <- cbind(rep(extended.row.idx, each = ext.col.len),
                      rep(extended.col.idx, times = ext.row.len))
  extended.scores <- SpatialExtremes::kriging(data = c(scores.mat),
                                      data.coord = cbind(rep(1:nr, nc),
                                                         rep(1:nc, each = nr)),
                                      krig.coord = krig.coord,
                                      cov.mod = "whitmat", grid = FALSE,
                                      sill = 1, range = 10,
                                      smooth = 0.8)$krig.est
  extended.scores <- matrix(extended.scores, nrow = ext.row.len,
                            ncol = ext.col.len, byrow = TRUE)

  # extended.scores = data.frame(extended.scores)
  extended.scores <- round(extended.scores, 3)

  row.dose <- as.numeric(rownames(scores.mat))
  col.dose <- as.numeric(colnames(scores.mat))

  extend.row.dose <- mapply(function(x, y){seq(from = x, to = y,
                                               length.out = len + 2)},
                            row.dose[-nr], row.dose[-1])
  extend.row.dose <- unique(round(c(extend.row.dose), 8))

  extend.col.dose <- mapply(function(x, y){seq(from = x, to = y,
                                               length.out = len + 2)},
                            col.dose[-nc], col.dose[-1])
  extend.col.dose <- unique(round(c(extend.col.dose), 8))

  rownames(extended.scores) <- extend.row.dose
  colnames(extended.scores) <- extend.col.dose

  return(extended.scores)
  # my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
  # heatmap.2(extended.scores, Rowv = F, Colv = F,
  #           dendrogram = 'none', trace = 'none',
  #           col = my_palette, density.info = "none")

  # Clean up
  gc()
}

# col.dose <- as.numeric(rownames(scores.mat))
# row.dose <- as.numeric(colnames(scores.mat))
#
# extend.row.dose <- mapply(function(x, y){seq(from = x, to = y,
#                                              length.out = len + 2)},
#                           row.dose[-nr], row.dose[-1])
# extend.row.dose <- unique(round(c(extend.row.dose), 7))
#
# extend.col.dose <- mapply(function(x, y){seq(from = x, to = y,
#                                              length.out = len + 2)},
#                           col.dose[-nc], col.dose[-1])
# extend.col.dose <- unique(round(c(extend.col.dose), 7))
#
# len.col <- length(extend.row.dose)
# len.row <- length(extend.col.dose)
# melt.mat <- reshape2::melt(scores.mat)
# test <- SpatialExtremes::kriging(data = melt.mat[, 3],
#                                 data.coord = as.matrix(melt.mat[, 1:2]),
#                                 krig.coord = cbind(rep(extend.col.dose,
#                                                        times = len.row),
#                                                    rep(extend.row.dose,
#                                                        each = len.col)),
#                                 cov.mod = "whitmat",
#                                 grid = FALSE,
#                                 sill = 1, range = 10,
#                                 smooth = 0.8)$krig.est
#
# test <- t(matrix(test, nrow = len.row,
#                ncol = len.col, byrow = TRUE))
#
# colnames(test) <- extend.row.dose
# rownames(test) <- extend.col.dose
