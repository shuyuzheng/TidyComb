#' Check the input template
#'
#' \code{CheckTemplate} function check the template inputted from users to
#' ensure they are acceptable to downstream functions.
#'
#' @details Function \code{CheckTemplate} checks the template file inputted by
#' users. If the template contains problems, an error will raise and a message
#' will given to user for checking their files.
#'
#' @param template a data frame. It is the input data from users for analysis
#'
#' @return Non
#' @author Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#' @export
#' @examples
#' data <- read.csv(system.file("template.csv", package = "TidyComb"),
#'                  stringsAsFactors = FALSE)
#' CheckTemplate(data)
#'
#' # Check some messy to data
#' # data_m <- data[, -2]
#' # CheckTemplate(data_m)
#'
CheckTemplate <- function(template){

  # Remove empty lines
  n_na <- (rowSums(is.na(template)) + rowSums(template == "", na.rm = TRUE))
  if (!all(n_na != ncol(template))) {
    message("Some empty rows were deleted from the uploaded file.")
    template <- template[n_na != ncol(template),]
  }

  # Column names
  missing_col <- setdiff(c("block_id", "drug_row", "drug_col", "inhibition",
                           "conc_r", "conc_c", "conc_r_unit", "conc_c_unit",
                           "cell_line_name", "drug_row", "drug_col"),
                         colnames(template))
  if (length(missing_col) != 0)
    stop('Column(s) "', paste(missing_col, collapse = '", "'),
         '" is(are) missing from uploaded data. Please check input file and',
         " re-upload the data.")

  # Check columns' data type
  dt <- sapply(template[, c("conc_c", "conc_r", "inhibition")], class)
  res <- names(dt)[!grepl("(integer|numeric)", dt)]
  if (length(res) != 0) {
    stop("The data type of column(s) '", paste(res, collapse = "', '"),
         "' must be 'numeric', 'integer'. Please check and re-upload ",
         "your file.")
  }

  # NA values is inavailable in columns block_id, Conc_r and Conc_c
  na <- apply(template[, c("block_id", "conc_r", "conc_c")], 2,
              function(x) sum(is.na(x)))
  if (sum(na) != 0) {
    m <- c("block_id","conc_r", "conc_c")[as.logical(na)]
    stop("There are missing values in column '", paste(m, collapse = "', '"),
         "'. Please check and re-upload file.")
  }

  # Duplicate drug_row, drug_col or cell line
  n <- template %>%
    dplyr::group_by(block_id) %>%
    dplyr::summarise(drug_row = dplyr::n_distinct(drug_row),
              drug_col = dplyr::n_distinct(drug_col),
              cell_line_name = dplyr::n_distinct(cell_line_name),
              conc_c_unit = dplyr::n_distinct(conc_c_unit),
              conc_r_unit = dplyr::n_distinct(conc_r_unit))
  m <- as.data.frame(n[, -1] > 1)
  m$block_id <-  as.character(n$block_id)
  rows <- m[apply(m[, -6], 1, sum) > 0, ]
  if (nrow(rows) != 0) {
    message <- apply(m, 1, function(x){
      comb <- colnames(m)[as.logical(c(x[1:5], FALSE))]
      if (length(comb) != 0){
        return(paste0("'", paste(comb, collapse = "', '"), "' of block ", x[6]))
      } else {
        return(NULL)
      }
    })
    if (is.list(message)) {
      message <- do.call(cbind, message)
    }
    stop("The columns 'drug_row', 'drug_col', 'cell_line_name', 'conc_c_unite',",
         " and 'conc_r_unite' can only have one value in one block. There are",
         " more than one values in:\n",
         paste0("colume(s) ", message, collapse = "; "), "\n",
         "Please check and re-upload the file.")
  }

  # Duplicate concentration combinations in one block
  block_dup <- duplicated(template[,c("conc_c", "conc_r", "block_id")])
  if (sum(block_dup) != 0) {
    stop("There are duplicated ('conc_r'-'conc_c') combinations in the block: ",
         paste(unique(template$block_id[block_dup]), collapse = ", "), ".\n",
         "Analysis function can not handle replicate data right now. Please ",
         "remove the duplicated rows and re-upload the data.")
  }

  gc()
  return(template)
}
