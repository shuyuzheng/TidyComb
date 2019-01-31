
#' calculate synergy scores in response table
#'
#' \code{GenerateScore()} takes input \code{template} table, calculates synergy
#' scores and generates the \em{response_with_score} table which is used for
#' downstream functions \code{\link{GenerateCurve()}},
#' \code{\link{GenerateSummary2()}}, and \code{\link{Generatesurface3()}}.
#'
#' There are four different synergy scores calculated by function
#' \code{GenerateScore}:
#' \itemize{
#'   \item \strong{Bliss}: Expected effect of two drugs acting indenpendently.
#'   \href{https://doi.org/10.1111/j.1744-7348.1939.tb06990.x}{(Bliss, 1939)}.
#'   \item \strong{Highest Single Agent(HSA)}: The maximal sigle drug effect.
#'   \href{http://pharmrev.aspetjournals.org/content/41/2/93.long}{(Berenbaum,
#'   1989)}.
#'   \item \strong{Loewe}: Expected effect of a drug combined with itself.
#'   \href{https://www.ncbi.nlm.nih.gov/pubmed/13081480}{(Loewe, 1953)}.
#'   \item \strong{Zero Interaction Potency(ZIP)}: Expected effect of two drugs
#'   that do not potentiate each other.
#'   \href{https://doi.org/10.1016/j.csbj.2015.09.001}{(Yadav.et.al., 2015)}
#' }
#'
#' @param response a data frame which must contains variables:
#' \itemize{
#'   \item \strong{block_id}: Numerical group ID, only positive. Identical value
#'   refers to screening results of one pair of compounds tested with all doses
#'   in one cancer cell line.
#'   \item \strong{conc_r}: Numerical value, only positive. Row compound
#'   concentration.
#'   \item \strong{conc_c}: Numerical value, only positive. Column compound
#'   concentration.
#'   \item \strong{response}: numerical value, positive or negative. Cell
#'   growth inhibition rate. It refers to the response of cancer cell lines to
#'   one drug combination.
#'   \item \strong{drug_row}: Names or identifiers of compounds which were added
#'   along rows.
#'   \item \strong{drug_col}: Names or identifiers of compounds which were added
#'   along columns.
#'   \item \strong{conc_r_unit}: Molar concentration unit of row compound with a
#'   relevant SI prefix.
#'   \item \strong{conc_c_unit}: Molar concentration unit of column compound
#'   with a relevant SI prefix.
#'   \item \strong{cell_line_name}: Name of cancer cell line in which compound
#'   combinations are tested.
#' }
#'
#' @return A data frame. Response table with synergy scores.
#'
#' @export
#'
#' @examples
#' response_with_score <- GenerateScore(system.file("extdata",
#'                                                  "template.csv",
#'                                                  package = "TidyComb"))
GenerateScore <- function(response, type = "inhibition"){
  if (!all(c("block_id", "drug_row", "drug_col", "response",
             "conc_r", "conc_c", "conc_r_unit", "conc_c_unit") %in%
           colnames(response))) {
    stop("The input data must contain the following columns: ",
         "block_id, drug_row, drug_col, response,","conc_r, conc_c, \n ",
         "conc_r_unit, conc_c_unit")
  }

  # Transform viability to inhibition
  if (type == "viability") {
    response$response <- 100 - response$response
  }

  # add row and column numbers
  response <- plyr::ddply(response, c("cell_line_name", "drug_row" ,
                                      "drug_col", "block_id"),
                          transform,
                          row = own_rank(conc_r),
                          col = own_rank(conc_c))

  scores <- list()
  m <- unique(response$block_id)

  # options(show.error.messages = F)
  for(i in m){

    data.tmp <- response[response$block_id == i,]

    set.seed(1) # add random noise - NB! the noise will be saved in the output
    data.tmp$response <- data.tmp$response + stats::rnorm(nrow(data.tmp),
                                                          0, 0.001)

    # response.mat <- ReshapeData2(data.tmp, data.type = "inhibition")
    # response.mat <- data.tmp2$dose.response.mats[[1]]

    response.mat <- reshape2::acast(data.tmp, conc_r ~ conc_c,
                                    value.var = "response",
                                    function(x) mean(x, na.rm = TRUE))
    # replicates are averaged

    # missing value imputation - NB! the imputation only stay within the
    # function
    # find the average of the neighboring cells in the matrix
    x <- array(c(rbind(response.mat[-1,], NA),
                 rbind(NA, response.mat[-nrow(response.mat), ]),
                 cbind(response.mat[,-1], NA),
                 cbind(NA, response.mat[, -ncol(response.mat)])),
               dim=c(nrow(response.mat), ncol(response.mat), 4))
    x.imp <- apply(x,c(1,2), function(x) mean(x, na.rm = TRUE))
    index.na <- is.na(response.mat)
    response.mat[index.na] <- x.imp[index.na]

    # one more round
    x <- array(c(rbind(response.mat[-1,], NA),
                 rbind(NA, response.mat[-nrow(response.mat),]),
                 cbind(response.mat[,-1], NA),
                 cbind(NA, response.mat[, -ncol(response.mat)])),
               dim=c(nrow(response.mat),ncol(response.mat), 4))
    x.imp <- apply(x, c(1,2), function(x) mean(x, na.rm=TRUE))
    index.na <- is.na(response.mat)
    response.mat[index.na] <- x.imp[index.na]

    # missing_index <- which(is.na(response.mat), arr.ind = T)
    # if(length(missing_index) !=0 ){
    #   for(j in 1:nrow(missing_index)){
    #     r <- missing_index[j, 1]
    #     c <- missing_index[j, 2]
    #
    #     tmp <- mean(c(response.mat[r + 1, c], response.mat[r - 1, c],
    #                  response.mat[r,c - 1], response.mat[r, c + 1]),
    #                na.rm = TRUE)
    #     if (is.na(tmp)) tmp <- 0 # if no neighbors are found
    #     response.mat[r, c] <- tmp
    #
    #   }
    # }

    # CalculateSynergy2 does not allow NA values
    hsa.tmp <- CalculateSynergy2(response.mat, method = "HSA", correction = TRUE,
                                 nan.handle = "L4", Emin = NA, Emax = NA)
    bliss.tmp <- CalculateSynergy2(response.mat, method = "BLISS",
                                   correction = TRUE, nan.handle = "L4",
                                   Emin = NA, Emax = NA)
    zip.tmp <- CalculateSynergy2(response.mat, method = "ZIP", correction = TRUE,
                                 nan.handle = "L4", Emin = NA, Emax = NA)

    loewe.tmp <- CalculateSynergy2(response.mat, method = "LOEWE",
                                   correction = TRUE, nan.handle = "L4",
                                   Emin = NA, Emax = NA)

    # data.tmp$synergy_zip <- apply(data.tmp[,c("row","col")], 1,
    #                               function(x) zip.tmp$scores[[1]][x[1],x[2]])
    # data.tmp$synergy_hsa <- apply(data.tmp[,c("row","col")], 1,
    #                               function(x) hsa.tmp$scores[[1]][x[1],x[2]])
    # data.tmp$synergy_bliss <- apply(data.tmp[,c("row","col")], 1,
    #                                 function(x) bliss.tmp$scores[[1]][x[1],x[2]])
    # data.tmp$synergy_loewe <- apply(data.tmp[,c("row","col")], 1,
    #                                 function(x) loewe.tmp$scores[[1]][x[1],x[2]])
    zip.tmp <- reshape2::melt(zip.tmp$scores)
    names(zip.tmp) <- c("conc_r", "conc_c", "synergy_zip")
    bliss.tmp <- reshape2::melt(bliss.tmp$scores)
    names(bliss.tmp) <- c("conc_r", "conc_c", "synergy_bliss")
    loewe.tmp <- reshape2::melt(loewe.tmp$scores)
    names(loewe.tmp) <- c("conc_r", "conc_c", "synergy_loewe")
    hsa.tmp <- reshape2::melt(hsa.tmp$scores)
    names(hsa.tmp) <- c("conc_r", "conc_c", "synergy_hsa")

    data.tmp <- full_join(data.tmp, zip.tmp, by = c("conc_r", "conc_c"))
    data.tmp <- full_join(data.tmp, hsa.tmp, by = c("conc_r", "conc_c"))
    data.tmp <- full_join(data.tmp, bliss.tmp, by = c("conc_r", "conc_c"))
    data.tmp <- full_join(data.tmp, loewe.tmp, by = c("conc_r", "conc_c"))

    scores[[i]] <-  data.tmp
  }
  # options(show.error.messages = TRUE)
  response_with_scores <- do.call(rbind, scores)
  return(response_with_scores)

}




# GenerateScore <- function(response, ...){
#   # add row and column numbers
#   response <- plyr::ddply(response, c("cell_line_name", "drug_row" ,
#                                       "drug_col", "block_id"),
#                           transform,
#                           row = own_rank(conc_r),
#                           col = own_rank(conc_c))
#   scores <- list()
#   m <- unique(response$block_id)
#   len <- length(m)
#
#   # options(show.error.messages = F)
#   for(i in 1:len){
#     message(round(i/len * 100), "%", "\r", appendLF = FALSE)
#     utils::flush.console()
#
#     index <- which(response$block_id == m[i])
#     data.tmp <- response[index,]
#
#     set.seed(1) # add random noise - NB! the noise will be saved in the output
#     data.tmp$response <- data.tmp$response + stats::rnorm(nrow(data.tmp),
#                                                           0, 0.001)
#
#     data.tmp2 <- ReshapeData2(data.tmp, data.type = "inhibition")
#     response.mat <- data.tmp2$dose.response.mats[[1]]
#
#     # missing value imputation - NB! the imputation only stay within the
#     # function
#     # find the average of the neighboring cells in the matrix
#     x <- array(c(rbind(response.mat[-1,], NA),
#                 rbind(NA, response.mat[-nrow(response.mat), ]),
#                 cbind(response.mat[,-1], NA),
#                 cbind(NA, response.mat[, -ncol(response.mat)])),
#               dim=c(nrow(response.mat), ncol(response.mat), 4))
#     x.imp <- apply(x,c(1,2), function(x) mean(x, na.rm = TRUE))
#     index.na <- is.na(response.mat)
#     response.mat[index.na] <- x.imp[index.na]
#
#     # one more round
#     x <- array(c(rbind(response.mat[-1,], NA),
#                 rbind(NA, response.mat[-nrow(response.mat),]),
#                 cbind(response.mat[,-1], NA),
#                 cbind(NA, response.mat[, -ncol(response.mat)])),
#               dim=c(nrow(response.mat),ncol(response.mat), 4))
#     x.imp <- apply(x, c(1,2), function(x) mean(x, na.rm=TRUE))
#     index.na <- is.na(response.mat)
#     response.mat[index.na] <- x.imp[index.na]
#
#     data.tmp2$dose.response.mats[[1]] <- response.mat
#
#     # missing_index <- which(is.na(response.mat), arr.ind = T)
#     # if(length(missing_index) !=0 ){
#     #   for(j in 1:nrow(missing_index)){
#     #     r <- missing_index[j, 1]
#     #     c <- missing_index[j, 2]
#     #
#     #     tmp <- mean(c(response.mat[r + 1, c], response.mat[r - 1, c],
#     #                  response.mat[r,c - 1], response.mat[r, c + 1]),
#     #                na.rm = TRUE)
#     #     if (is.na(tmp)) tmp <- 0 # if no neighbors are found
#     #     response.mat[r, c] <- tmp
#     #
#     #   }
#     # }
#
#     # CalculateSynergy2 does not allow NA values
#     hsa.tmp <- CalculateSynergy2(data.tmp2, method = "HSA", correction = TRUE,
#                                  nan.handle = "L4", Emin = NA, Emax = NA)
#     bliss.tmp <- CalculateSynergy2(data.tmp2, method = "BLISS",
#                                    correction = TRUE, nan.handle = "L4",
#                                    Emin = NA, Emax = NA)
#     zip.tmp <- CalculateSynergy2(data.tmp2, method = "ZIP", correction = TRUE,
#                                 nan.handle = "L4", Emin = NA, Emax = NA)
#     loewe.tmp <- CalculateSynergy2(data.tmp2, method = "LOEWE",
#                                   correction = TRUE, nan.handle = "L4",
#                                   Emin = NA, Emax = NA)
#
#     data.tmp$synergy_zip <- apply(data.tmp[,c("row","col")], 1,
#                                   function(x) zip.tmp$scores[[1]][x[1],x[2]])
#     data.tmp$synergy_hsa <- apply(data.tmp[,c("row","col")], 1,
#                                  function(x) hsa.tmp$scores[[1]][x[1],x[2]])
#     data.tmp$synergy_bliss <- apply(data.tmp[,c("row","col")], 1,
#                                    function(x) bliss.tmp$scores[[1]][x[1],x[2]])
#     data.tmp$synergy_loewe <- apply(data.tmp[,c("row","col")], 1,
#                                    function(x) loewe.tmp$scores[[1]][x[1],x[2]])
#
#     scores[[i]] <- data.tmp
#   }
#   # options(show.error.messages = TRUE)
#   response_with_scores <- do.call(rbind, scores)
#   return(response_with_scores)
#
# }


