################################################################################
# Copyright Shuyu Zheng and Jing Tang - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Shuyu Zheng <shuyu.zheng@helsinki.fi>, November 2020
################################################################################

# TidyComb
# Functions for calculating drug combination response or synergy surface
#
# Functions on this page:
# smoothing

#' Smoothing the dose response or synergy score surface
#'
#' This function will add the data points into the existing score matrix to make
#' the surface in 3D visualization smoothly. Kriging algotithom is used for
#' imputation.
#'
#' @param scores.mat matrix, the matix for response valuses or synergy scores.
#' @param len integer, the number of data points will be added between two
#' neighbor dosages in original matrix.
#'
#' @return The matrix with addational data points.
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
smoothing <-  function (scores.mat, len = 2) {

  options(scipen = 999)
  nr <- nrow(scores.mat)
  nc <- ncol(scores.mat)

  # missing value imputation
  while (sum(is.na(scores.mat))) {
    scores.mat <- synergyfinder::ImputeNA(scores.mat)
  }
  ext.row.len <- (nr - 1) * (len + 2) - (nr - 2)
  ext.col.len <- (nc - 1) * (len + 2) - (nc - 2)

  extended.row.idx <- seq(1, nr, length = ext.row.len)
  extended.col.idx <- seq(1, nc, length = ext.col.len)

  krig.coord <- cbind(rep(extended.row.idx, each = ext.col.len),
                      rep(extended.col.idx, times = ext.row.len))
  extended.scores <- synergyfinder::kriging(data = c(scores.mat),
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
  extend.row.dose <- c(extend.row.dose[1, 1], extend.row.dose[-1, ])

  extend.col.dose <- mapply(function(x, y){seq(from = x, to = y,
                                               length.out = len + 2)},
                            col.dose[-nc], col.dose[-1])
  extend.col.dose <- c(extend.col.dose[1, 1], extend.col.dose[-1, ])

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
