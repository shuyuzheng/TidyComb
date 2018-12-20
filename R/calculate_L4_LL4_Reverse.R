# Reverse function of L.4
L4Reverse <- function(b, c, d, e, y) {
  x <- (log((d - c) / (y - c) - 1)) / b + e
  return(x)
}

# Reverse function of LL.4
LL4Reverse <- function(b, c, d, e, y) {
  x <- e * (((d - c) / (y - c) - 1) ^ (1 / b))
  return(x)
}