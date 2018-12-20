FittingSingleDrug2 <- function (response.mat,
                                fixed = c(NA, NA, NA, NA),
                                nan.handle = c("L4")) {
  r.num <- nrow(response.mat)
  c.num <- ncol(response.mat)
  drug.col <- cbind(as.numeric(colnames(response.mat)),
                    response.mat[1, 1:c.num])
  drug.col <- data.frame(drug.col)
  colnames(drug.col) <- c("conc", "effect")
  drug.col$conc <- as.numeric(drug.col$conc)
  # nonzero concentrations to take the log
  drug.col$conc[which(drug.col$conc==0)] <- 10^-10
  drug.col$effect <- as.numeric(drug.col$effect)

  if(nrow(drug.col)!= 1) {
    if (stats::var(drug.col$effect) == 0) {
      drug.col$effect[nrow(drug.col)] <- drug.col$effect[nrow(drug.col)] +
        10^-10
    }
  }

  nan.handle <- match.arg(nan.handle)
  drug.col.model <- tryCatch({
    drc::drm(effect ~ conc, data = drug.col, fct = drc::LL.4(fixed = fixed),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  }, warning = function(w) {
    drc::drm(effect ~ log(conc), data = drug.col, fct = drc::L.4(fixed = fixed),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  }, error = function(e) {
    drc::drm(effect ~ log(conc), data = drug.col, fct = drc::L.4(fixed = fixed),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  })
  drug.col.fitted <- suppressWarnings(stats::fitted(drug.col.model))


  drug.row <- cbind(as.numeric(rownames(response.mat)),
                    response.mat[1:r.num, 1])
  drug.row <- data.frame(drug.row)
  colnames(drug.row) <- c("conc", "effect")

  drug.row$conc = as.numeric(drug.row$conc)
  # nonzero concentrations to take the log
  drug.row$conc[which(drug.row$conc==0)] <- 10^-10
  drug.row$effect <- as.numeric(drug.row$effect)

  if (nrow(drug.row)!=1){
    if (stats::var(drug.row$effect) == 0) {
     drug.row$effect[nrow(drug.row)] <- drug.row$effect[nrow(drug.row)] +
        10^-10
    }
  }
  drug.row.model <- tryCatch({
    drc::drm(effect ~ conc, data = drug.row,
             fct = drc::LL.4(fixed = fixed),
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  }, warning = function(w) {
    drc::drm(effect ~ log(conc),
             data = drug.row,
             fct = drc::L.4(fixed = fixed),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  }, error = function(e) {
     drc::drm(effect ~ log(conc),
              data = drug.row,
              fct = drc::L.4(fixed = fixed),
              na.action = stats::na.omit,
              control = drc::drmc(errorm = FALSE, noMessage = TRUE))

  })
  drug.row.fitted <- suppressWarnings(stats::fitted(drug.row.model))

  res <- list(drug.row.fitted = drug.row.fitted,
             drug.row.model = drug.row.model,
             drug.col.model = drug.col.model,
             drug.col.fitted = drug.col.fitted)
  return(res)
}
