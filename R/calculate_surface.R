# TidyComb
# Functions for Calculating drug combination response or sysurface
# Copyright Shuyu Zheng
#
# Functions on this page:
#
smoothing <-  function (scores.mat, len = 2) {

  options(scipen = 999)
  nr <- nrow(scores.mat)
  nc <- ncol(scores.mat)

  # missing value imputation
  score.mat <- ImputeNear(scores.mat)

  ext.row.len <- (nr - 1) * (len + 2) - (nr - 2)
  ext.col.len <- (nc - 1) * (len + 2) - (nc - 2)

  extended.row.idx <- seq(1, nr, length = ext.row.len)
  extended.col.idx <- seq(1, nc, length = ext.col.len)

  krig.coord <- cbind(rep(extended.row.idx, each = ext.col.len),
                      rep(extended.col.idx, times = ext.row.len))
  extended.scores <- SpatialExtremes::kriging(data = scores.mat,
                                      data.coord = cbind(rep(1:nr, times = nc),
                                                         rep(1:nc, each = nr)),
                                      krig.coord = krig.coord,
                                      cov.mod = "whitmat", grid = FALSE,
                                      sill = 1, range = 10,
                                      smooth = 0.8)$krig.est
  extended.scores <- matrix(extended.scores, nrow = ext.row.len,
                            ncol = ext.col.len, byrow = TRUE)

  # extended.scores = data.frame(extended.scores)
  extended.scores <- round(extended.scores, 3)

  col.dose <- as.numeric(colnames(scores.mat))
  row.dose <- as.numeric(rownames(scores.mat))

  row.index <- 1:nr
  col.index <- 1:nc

  row.tick <- which(extended.row.idx %in% row.index)
  col.tick <- which(extended.col.idx %in% col.index)

  rownames(extended.scores)[row.tick] <- row.dose
  colnames(extended.scores)[col.tick] <- col.dose

  list1 <- lapply(diff(row.dose)/diff(row.tick),
                  function(x) cumsum(rep(x, each = len)))
  list2 <- as.list(row.dose[-nr])
  rownames(extended.scores)[-row.tick] <- as.character(round(unlist(mapply("+",
                                                                           list1, list2, SIMPLIFY = FALSE)), 7))

  list1 <- lapply(diff(col.dose)/diff(col.tick),
                  function(x) cumsum(rep(x, each = len)))
  list2 <-  as.list(col.dose[-nc])
  colnames(extended.scores)[-col.tick] <- as.character(round(unlist(mapply("+",
                                                                           list1, list2, SIMPLIFY = FALSE)), 7))

  return(extended.scores)
  # my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
  # heatmap.2(extended.scores, Rowv = F, Colv = F,
  #           dendrogram = 'none', trace = 'none',
  #           col = my_palette, density.info = "none")

}