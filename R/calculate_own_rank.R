#' Rank the vectors
#'
#' @param x a vector need to be ranked
#'
#' @return a vector contains the rank of original vector
#' @export
own_rank <- function(x) {
  x_unique <- unique(x)
  x_ranks <- rank(x_unique)
  x_ranks[match(x, x_unique)]
}

# own_log, natural base
# calculate log(1+10^(b*(c-x))) to be used in scoreCurve function
own_log <-  function(b, c, x) {
  arg = 1 + 10 ^ (b * (c - x))
  if(is.infinite(arg) == TRUE) {
    res <- b * (c - x) * log(10)
  } else {
    res <- log(arg)
  }
  return(res)
}

# own_log2, natural base
# calculate log(1+exp(x)) to be used in scoreCurve.L4 function
own_log2 <- function(x) {
  arg <- 1 + exp(x)
  if(is.infinite(arg) == TRUE) {
    res <-  x
  } else {
    res = log(arg)
  }
  return(res)
}
