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
