
#' Imputation with nearest values
#'
#' \code{ImputeNear} do missing value inputation by averaging 4 values in
#' nearest cells (top, bottom, left, right) to the cell comtains NA value.
#'
#' @param response.mat A matrix which has missing value.
#'
#' @param times a integer which refer to the time of runing imputation.
#'
#' @return A matrix which is same as input matrix except the NA cells are filled
#' with numbers.
#'
#' @export
ImputeNear <- function(response.mat, times = 1) {
  for (i in seq(to = times)) {
  x <- array(c(rbind(response.mat[-1,], NA),
               rbind(NA, response.mat[-nrow(response.mat), ]),
               cbind(response.mat[,-1], NA),
               cbind(NA, response.mat[, -ncol(response.mat)])),
             dim=c(nrow(response.mat), ncol(response.mat), 4))
  x.imp <- apply(x,c(1,2), class)# function(x) mean(x, na.rm = TRUE))
  index.na <- is.na(response.mat)
  response.mat[index.na] <- x.imp[index.na]
  }
  return(response.mat)
}
