# Function to generate the response table with synergy scores

# response input must contain the following columns
# block_id	conc_r	conc_c	response	drug_row	drug_col	conc_r_unit
# conc_c_unit	cell_line_name

GenerateScore <- function(response){
  # add row and column numbers
  response <- plyr::ddply(response, c("cell_line_name", "drug_row" ,
                                      "drug_col", "block_id"),
                          transform,
                          row = own_rank(conc_r),
                          col = own_rank(conc_c))
  scores <- list()
  m <- unique(response$block_id)
  len <- length(m)

  # options(show.error.messages = F)
  for(i in 1:len){
    # cat('\r', i)
    print(m[i])
    index <- which(response$block_id == m[i])
    data.tmp <- response[index,]

    set.seed(1) # add random noise - NB! the noise will be saved in the output
    data.tmp$response <- data.tmp$response + stats::rnorm(nrow(data.tmp),
                                                          0, 0.001)

    data.tmp2 <- ReshapeData2(data.tmp, data.type = "inhibition")
    response.mat <- data.tmp2$dose.response.mats[[1]]

    # missing value imputation - NB! the imputation only stay within the function
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

    data.tmp2$dose.response.mats[[1]] <- response.mat

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
    hsa.tmp <- CalculateSynergy2(data.tmp2, method = "HSA", correction = TRUE,
                                 nan.handle = "L4", Emin = NA, Emax = NA)
    bliss.tmp <- CalculateSynergy2(data.tmp2, method = "BLISS",
                                   correction = TRUE, nan.handle = "L4",
                                   Emin = NA, Emax = NA)
    zip.tmp <- CalculateSynergy2(data.tmp2, method = "ZIP", correction = TRUE,
                                nan.handle = "L4", Emin = NA, Emax = NA)
    loewe.tmp <- CalculateSynergy2(data.tmp2, method = "LOEWE",
                                  correction = TRUE, nan.handle = "L4",
                                  Emin = NA, Emax = NA)

    data.tmp$synergy_zip <- apply(data.tmp[,c("row","col")], 1,
                                  function(x) zip.tmp$scores[[1]][x[1],x[2]])
    data.tmp$synergy_hsa <- apply(data.tmp[,c("row","col")], 1,
                                 function(x) hsa.tmp$scores[[1]][x[1],x[2]])
    data.tmp$synergy_bliss <- apply(data.tmp[,c("row","col")], 1,
                                   function(x) bliss.tmp$scores[[1]][x[1],x[2]])
    data.tmp$synergy_loewe <- apply(data.tmp[,c("row","col")], 1,
                                   function(x) loewe.tmp$scores[[1]][x[1],x[2]])

    scores[[i]] <- data.tmp
  }
  # options(show.error.messages = TRUE)
  response_with_scores <- do.call(rbind, scores)
  return(response_with_scores)

}


