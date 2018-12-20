d2.fun <- function(col_conc, row_conc, drug.row.model) {
  drug.row.par <- stats::coef(drug.row.model)
  # LL.4, conc must be raw
  if(length(grep("LL.4", drug.row.model$call$fct)) > 0 ) {
    conc = col_conc + row_conc
    (drug.row.par[3] + drug.row.par[2] *
        (conc / drug.row.par[4]) ^ drug.row.par[1]) /
      (1 + (conc / drug.row.par[4]) ^ drug.row.par[1])
  } else {#L.4, conc must be logscaled, ie. log(conc)
    conc = log(col_conc + row_conc)
    (drug.row.par[2] + (drug.row.par[3] - drug.row.par[2]) /
      (1 + exp(drug.row.par[1] * (conc - drug.row.par[4]))))
  }
}