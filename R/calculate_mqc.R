#' Calculate the mQC descriptors from a combination response matrix
#'
#' @param mat A numeric, square matrix (required)
#' @param reverse A boolean flag indicating the type of assay readout. A typical
#'   combination screening readout such as cell titer glow should have negative
#'   control (DMSO) position at [N,N] and have signal equal to 100. However,
#'   when dealing with other readouts that have negative control signal equal to
#'   0 (e.g., Caspase-Glo), set reverse = \code{TRUE} (default = \code{FALSE})
#' @param opt A vector that controls the contents in the output list. available
#'   option: dmso, sa, dsa, moran, smoothness, mono. (default = c('dmso', 'sa',
#'   'moran', 'smoothness', 'mono'))
#' @param permutation A integer indicating the number of permutations in the
#'   moothness test (default = 10000)
#' @param coeff A numeric vector used in dsa calculation (deprecated,
#'   default = c(0.8811, 8.8209))
#' @param fullmatrix If \code{TRUE}, Moran I test will be performed on the full
#'   matrix; if \code{FALSE}, Moran I test will be performed on the combination
#'   submatrix. (default = "TRUE")
#' @param ncpu A numeric integer indicating the number of processors.
#'   (default = 1)
#' @param na.action A string indicating how to treat NA in the matrix. "bad":
#'   return a vector leading to "Bad" prediction; "noactivity": substitute 100
#'   or 0 if reverse=\code{TRUE} (default = "bad")
#' @return A vector of feature values
#' @author Lu Chen \email{chenl22@@mail.nih.gov}
#' @references Chen, L. et al., "mQC: A Heuristic Quality-Control Metric for
#'   High-Throughput Drug Combination Screening", Sci. Rep., 2016, 6:37741
#' @export

resp.qc.calc <- function(mat, reverse=FALSE, opt=c('dmso', 'sa', 'moran',
                                                   'smoothness', 'mono'),
                         permutation=10000, coeff=c(0.8811, 8.8209),
                         fullmatrix=TRUE, ncpu = 1, na.action="bad") {
  status <- require(mgcv) && require(ape) &&
            ((ncpu > 1 && require(parallel)) || ncpu == 1)
  if (!status) stop("mQC requires mgcv and ape package. \n",
                    "The multiprocessing version requires parallel package")
  n <- dim(as.matrix(mat))[1]
  r <- as.matrix(mat)

  ##########################
  # 0.initialization (bad) #
  ##########################
  dmso.v <- 0
  sa.min <- 0
  sa.max <- 0
  sa.matrix <- 0
  dsa.v <- 1e3    # deprecated feature
  moran.p <- 1
  smoothness.p <- 1
  mono.v <- 0

  ## intitialize return vector
  rst <- c()

  ## parse options
  opt.st <- c('dmso', 'sa', 'dsa.use', 'moran', 'smoothness', 'mono')
  sel <- opt.st %in% opt
  dmso.use <- sel[1]
  sa.use <- sel[2]
  dsa.use <- sel[3]
  moran.use <- sel[4]
  smoothness.use <- sel[5]
  mono.use <- sel[6]

  ## NA actions
  if (!na.action %in% c("bad", "noactivity")) stop("na.action must be 'bad' ",
                                                   "or 'noactivity'")
  if (na.action == "bad" && any(is.na(r))) {
    if (dmso.use) rst <- c(rst, dmso.v=dmso.v)
    if (sa.use) rst <- c(rst, sa.min=sa.min, sa.max=sa.max, sa.matrix=sa.matrix)
    if (dsa.use) rst <- c(rst, dsa.v=dsa.v)
    if (moran.use) rst <- c(rst, moran.p=moran.p)
    if (smoothness.use) rst <- c(rst, smoothness.p=smoothness.p)
    if (mono.use) rst <- c(rst, mono.v=mono.v)
    return(rst)
  }
  if (na.action == "noactivity") r[is.na(r)] <- ifelse(reverse, 0, 100)

  ## do QC on bounded data
  r[r<0] <- 0
  r[r>100] <- 100

  ####################
  # 1. dmso response #
  ####################
  ## check if dmso response is good
  ## may change to plate quality metric, e.g., z-score in the future
  if (dmso.use) {
    dmso.v <- ifelse(reverse, 100 - r[n,n], r[n,n])
    rst <- c(rst, dmso.v = dmso.v)
  }

  ##################################
  # 2. relative standard deviation #
  ##################################
  ## Use relative SD http://en.wikipedia.org/wiki/Relative_standard_deviation
  rsd <- function(x)  {
    ret <- sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE) * 100
    if (is.na(ret) || is.nan(ret)) return(0)
    return(ret)
  }

  if (sa.use) {
    sa.min <- rsd(r[,n])
    sa.max <- rsd(r[n,])
    sa.matrix <- rsd(as.numeric(r[-n,-n]))
    # swap sa.min and sa.max
    if (sa.min > sa.max) {
      tmp <- sa.min
      sa.min <- sa.max
      sa.max <- tmp
    }
    rst <- c(rst, sa.min = sa.min, sa.max = sa.max, sa.matrix = sa.matrix)
  }

  #####################################
  # 3. delta rsd based on bliss model #
  #####################################
  # res.matrix = res.A * res.B / 100
  # coeff should be provided
  # deprecated
  if (sa.use && dsa.use) {
    dsa.v <- sqrt(sa.min^2 + sa.max^2) - (coeff[1] * sa.matrix + coeff[2])
    rst <- c(rst, dsa.v = dsa.v)
  }

  ##############################
  # 4. spatial autocorrelation #
  ##############################
  ## http://www.ats.ucla.edu/stat/r/faq/morans_i.htm
  ## check spatial autocorrelation using moran's I test
  moran <- function(r, fullmatrix = fullmatrix) {
    n <- dim(r)[1]
    if (fullmatrix) {
      sr <- as.numeric(r)
    } else {
      sr <- as.numeric(r[-n, -n])
      n <- n - 1
    }
    if (length(unique(sr)) == 1) {
      mp <- 1
    } else {
      sr <- cbind(expand.grid(1:n, 1:n), sr)
      srd <- 1/as.matrix(dist(sr[,1:2]))
      diag(srd) <- 0
      mtest <- ape::Moran.I(sr[,3], srd, na.rm=TRUE)
      return(mtest$p.value)
    }
  }

  ## calculate the moran's I p-value
  if (moran.use) {
    moran.p <- moran(r, fullmatrix=fullmatrix)
    rst <- c(rst, moran.p = moran.p)
  }

  ##########################################
  # 5.smoothness of the activity landscape #
  ##########################################
  ## calculate the rmsd between the smoothed grids values and original ones
  ## using generalized additive models (gam)
  ## metrics to measure the smoothness of the activity landscape
  ## Note: worker for smoothness.test function
  ## 10/26/2015: will throw error if the matrix is too small (w/o NA omit process)
  smoothness.calc <- function(m) {
    N <- nrow(m)
    status <- try(
      mod <- mgcv::gam(z ~ te(x, y), na.action=na.omit,
                       data = data.frame(x = rep(seq_len(N), each = N),
                                         y = rep(seq_len(N), times = N),
                                         z = c(m))),
      silent=TRUE
    )
    if (inherits(status, 'try-error')) {
      return(NA)
    }
    ## some NA may be omit by gam optimizor
    rmsd <- sqrt(sum((mod$residuals)^2)/ length(mod$residuals))
    return(rmsd)
  }

  ## perform statistical testing and return significance of smoothness
  ## permutation numerically, require(parallel) package for parallelization
  smoothness.test <- function(m, permutation=permutation, ncpu=ncpu,
                              na.action=na.omit) {
    n <- nrow(m)
    ## for flat landscape, return p-value = 0
    if (min(na.action(m)) == max(na.action(m))) return(0)
    ## if smoothness metric is non-calculateable, return p-value = 1
    rmsd <- smoothness.calc(m)
    if (is.na(rmsd)) return(1)
    ## shuffle matrix numerically
    shuffle_matrix <- function(x, m) {
      m2 <- sample(m)
      n <- nrow(m)
      dim(m2) <- c(n,n)
      smth <- smoothness.calc(m2)
      return(smth)
    }
    ## start statistical testing - generate the reference population
    if (ncpu == 1) {
      pop <- lapply(seq(permutation), shuffle_matrix, m = m)
    } else {
      pop <- do.call(c, parallel::mclapply(seq(permutation), shuffle_matrix,
                                           m=m, mc.cores=ncpu))
    }
    ## calculate the p-value based on numerically calculated cdf(pop)
    pop <- na.omit(unlist(pop))
    fprob <- ecdf(pop)
    smooth.pvalue <- fprob(rmsd)
    return(smooth.pvalue)
  }

  ## calculate smoothness p-value
  if (smoothness.use) {
    smoothness.p <- smoothness.test(r, ncpu=ncpu, permutation=permutation)
    rst <- c(rst, smoothness.p = smoothness.p)
  }

  #############################################
  # 6. monotinicity of the activity landscape #
  #############################################
  # for each response, calculate the number of items in the top-left and
  # bottem-right submatices which satisfy the monotonicity criteria, then
  # aggregate and normalize
  mono <- function(m) {
    rst <- c()
    n <- dim(m)[1]
    for (i in 1:n) {
      for (j in 1:n) {
        if (is.na(m[i,j]) || is.nan(m[i,j])) {
          val <- NA
        } else {
          if (!reverse) {
            lefttop <- sum(m[1:i, 1:j] <= m[i,j], na.rm=TRUE)
            rightbottom <- sum(m[i:n, j:n] >= m[i,j], na.rm=TRUE)
          } else {
            lefttop <- sum(m[1:i, 1:j] >= m[i,j], na.rm=TRUE)
            rightbottom <- sum(m[i:n, j:n] <= m[i,j], na.rm=TRUE)
          }
          val <- (lefttop + rightbottom - 2) / (i*j + (n-i+1)*(n-j+1) - 2)
        }
        rst <- c(rst, val)
      }
    }
    return(mean(rst, na.rm=TRUE))
  }

  ## calculate the monotinicty value
  if (mono.use) {
    mono.v <- mono(r)
    rst <- c(rst, mono.v = mono.v)
  }

  ##############
  # 7. output  #
  ##############
  return (rst)
}

#' Predict the mQC label and confidence from a combination response matrix or
#' pre-computed descriptors this function requires an adaboost model object and
#' a confidence model to make prediction by default, the model object is in
#' 'qcModel.rda', and the confidence model is in 'confidenceModel.rda'
#'
#' @param m A numeric, square response matrix. A list of matrices or blocks are
#'   also valid (either m or df is required)
#' @param df A numeric dataframe of pre-computed mQC descriptors. If not
#'   \code{NULL}, the mQC label will be predicted using the features in this
#'   data frame (either m or df is required default = NULL)
#' @param qcModel A adaboost model object for class prediction. If \code{NULL},
#'   the default qcModel will be loaded (default = NULL)
#' @param confidenceModel A dataframe for confidence calculation. If
#'   \code{NULL}, the default confidenceModel will be loaded (default = NULL)
#' @param ... Other arguments that pass to \code{resp.qc.calc} function
#' @return A data frame including the probability of being Good, Medium and Bad,
#'   the class and its confidence
#' @seealso \code{\link{resp.qc.calc}}
#' @author Lu Chen \email{chenl22@@mail.nih.gov}
#' @examples
#' # a good block
#' block1 <- matrix(c(5.43,6.57,6.08,7.29,6.75,7.20,9.84,7.75,7.43,6.12,2.59,
#'   2.91,2.47,2.96,2.59,2.28,3.31,4.24,4.29,3.16,3.50,2.20,3.21,2.57,2.82,2.17,
#'   3.87,3.53,3.22,5.54,4.09,4.00,2.47,3.11,5.12,5.22,4.47,5.06,3.80,6.60,5.57,
#'   6.26,4.42,5.49,8.20,8.34,10.10,9.84,10.69,10.57,31.80,34.99,46.76,57.36,
#'   58.36,74.10,85.89,82.35,79.52,94.14,89.38,78.18,91.26,100.00,89.28,92.86,
#'   86.12,83.62,100.00,91.78,81.88,85.73,92.95,100.00,83.82,100.00,100.00,
#'   100.00,100.00,94.90,77.19,100.00,84.36,100.00,96.38,100.00,100.00,100.00,
#'   100.00,87.47,72.38,89.91,100.00,100.00,100.00,100.00,90.17,100.00,92.43,
#'   98.95), nrow=10)
#' # a random block by shuffling block1
#' block2 <- matrix(sample(c(block1)), nrow=10)
#' # calculate block.qc for block1 and block2 using \code{block.qc} function
#' qc1 <- block.qc(m = list(block1, block2), permutation = 500, ncpu = 2)
#' # calculate mQC descriptor using \code{resp.qc.calc} first, and then
#' # calculate mQC with the descriptors
#' descriptor <- do.call(rbind, lapply(list(block1, block2), resp.qc.calc,
#'                       permutation = 500))
#' qc2 <- block.qc(df = descriptor)
#' # compare qc1 and qc2 (should be identical)
#' print(qc1)
#' print(qc2)
#' @export

block.qc <- function(m=NULL, df=NULL, qcModel=NULL, confidenceModel=NULL, ...) {
  #######################
  ## utility functions ##
  #######################
  # validate the data frame before proceeding
  validateDescriptorDf <- function(df) {
    df <- try(df[, c("dmso.v", "sa.min", "sa.max", "moran.p", "sa.matrix",
                     "smoothness.p", "mono.v")], silent = T)
    if (inherits(df, 'try-error')) {
      stop('Incorrect attribute names in the descriptor dataframe.')
    }
    return(df)
  }

  # extrapolation function
  extrapolate <- function(x1, x2, y1, y2, x) {
    if (x1 != x2) {
      y <- x* (y2-y1)/(x2-x1) + (x2*y1-x1*y2)/(x2-x1)
    } else {
      y <- ifelse((x1 == x2 && y1 == y2 && x1 == x), y1, NA)
    }
    return(y)
  }

  # calculate the confidence
  confidence.calc <- function(probs, model) {
    # make sure the model level are sorted
    model <- model[order(model$level), ]
    confidence <- apply(probs, 1, function(x) {
      val <- sd(x)
      # find the corresponding level
      foo <- abs(model$level - val)
      idx <- which(foo==min(foo))

      # calculate the corresponding confidence
      if (length(idx) > 1) {
        return(mean(model[idx, ]$confidence))
      }
      if (val == model[idx, ]$level) {
        return(model[idx, ]$confidence)
      }
      if (val > model[idx, ]$level) {
        if (idx == nrow(model)) {
          return(model[idx, ]$confidence)
        } else {
          return(extrapolate(x1 = model[idx, ]$level,
                             x2 = model[idx+1, ]$level,
                             y1 = model[idx, ]$confidence,
                             y2 = model[idx+1, ]$confidence,
                             x = val))
        }
      }
      if (val < model[idx, ]$level) {
        if (idx == 1) {
          return(model[idx, ]$confidence)
        } else {
          return(extrapolate(x1 = model[idx, ]$level,
                             x2 = model[idx-1, ]$level,
                             y1 = model[idx, ]$confidence,
                             y2 = model[idx-1, ]$confidence,
                             x = val))
        }
      }
    })
    return(confidence)
  }

  ## main function
  status <- require(adabag)
  if (!status) stop("mQC prediction requires adabag package")

  if (!is.null(df)) {
    descriptor <- validateDescriptorDf(df)
  } else {
    if (!is.list(m)) {
      m <- list(m)
    }
    descriptor <- do.call(rbind, lapply(m, resp.qc.calc, ...))
  }
  if (is.null(qcModel) || is.na(qcModel)) load('qcModel.rda')
  if (is.null(confidenceModel) || is.na(confidenceModel)) load('confidenceModel.rda')

  descriptor <- as.data.frame(descriptor)
  pred <- adabag::predict.boosting(qcModel, newdata = descriptor)

  # return a list
  invisible(data.frame(prob = pred$prob,
                       class = pred$class,
                       confidence = confidence.calc(pred$prob, confidenceModel)))
}